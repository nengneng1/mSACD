% =========================================================================
% step2_TMI_rlsu.m — PCA去噪 + TMI 通道拆分 (Richardson-Lucy 光谱解混)
% =========================================================================
% 输入 (workspace): stack1_for_TMI, imgstack, frame, paddingfactor,
%                   FILE_OUT_DIR, fname, confname
% 输出 (workspace): amp_raw, comp_raw, TMIch1, TMIch2, H, W
% 输出 (tif):       {FILE_OUT_DIR}/TMI_rlsu/
%
% 解混原理 (Kumar et al., Nat. Photonics 2025, RLSU):
%   线性混合模型:  decay(t) ≈ c1*ch1(t) + c2*ch2(t)
%   写成矩阵:      D = M * C,   M=[ch1, ch2] (F×2),  C=[c1; c2] (2×N)
%   Richardson-Lucy 迭代 (乘法更新, 泊松似然下推导):
%       C <- C .* ( M' * ( D ./ (M*C) ) ) ./ ( M' * 1 )
%   - 初值全 1, 乘法更新 → 估计值恒为非负 (无需截断/约束)
%   - 天然处理泊松散粒噪声, 低 SNR 数据更稳健
%   - 全向量化 (按像素分块), 无逐像素 fmincon, 比 fmincon 快几个数量级
% =========================================================================

%% ====== 路径 ======
STEP_NAME = 'TMI_rlsu';
UNMIX_TAG = 'rlsu';      % 方法标签 → 供 step2b/step3 区分输出目录
addpath(genpath('F'));

%% ====== 参数 ======
bg_thresh  = 50/65535;
pca_on     = 1;          % PCA 去噪开关 (0=关闭, 用于对比)
num_iters  = 1000;        % RL 迭代次数 (2 组分收敛快, 50~200 即可)
chunk_px   = 50000;      % 每块处理的像素数 (控制内存; 越大越快越占内存)

% --- Ratio-map TV 正则 (ROF, Chambolle; 与 phasor 版一致) ---
%   只对 c1 比例图做一次 TV, c2 = 1 - c1 (sum-to-one 下两通道梯度相反, TV 等价)
%   作用对象是"解混后的标量比例场"(分片常数+锐边界), 保边界、压平内部噪声
tv_on     = 1;           % 1=对 c1 比例图做 TV 正则, 0=关闭
tv_weight = 0.02;         % 正则强度 λ (越大越平滑; 0.02~0.1 常用)
tv_iter   = 100;         % Chambolle 迭代次数 (2D ROF 收敛快)
% save_iters = [1 5 20 50 100 500 1000];   % 额外保存这些迭代步的比例图, 便于观察收敛
save_iters = [1 5 20 50 75 100 200 300 1000];   % 额外保存这些迭代步的比例图, 便于观察收敛
%% ====== 输出目录 ======
out_dir = fullfile(FILE_OUT_DIR, STEP_NAME);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ====== 衰减曲线 ======
load([confname, '.mat'], 'ch1', 'ch2');
ch1 = real(ch1(:));  ch2 = real(ch2(:));   % 确保列向量、实数

%% ====== 参考曲线 (无 gamma, 保持线性混合假设) ======
stack1_g = stack1_for_TMI; ch1_g = ch1; ch2_g = ch2;

%% ====== 参考曲线可分性诊断 ======
R_ref = [ch1_g(1:frame), ch2_g(1:frame)];
c_ref = corrcoef(R_ref);
fprintf('  Reference: corr(ch1,ch2)=%.4f, cond([ch1,ch2])=%.2f\n', ...
    c_ref(1,2), cond(R_ref));

%% ====== PCA 去噪 (与 phasor/fmincon 版本一致; pca_on=0 时跳过) ======
stack1_norm = double(stack1_g ./ max(stack1_g(:)));
if pca_on == 1
    filt_no_PC  = 3;
    NoiseCorrection = sqrt(mean(stack1_norm, [1, 2]) + eps);
    temp_norm = stack1_norm ./ NoiseCorrection;
    temp_norm(temp_norm < 0) = 0;
    vector_conv = reshape(temp_norm, [size(temp_norm,1)*size(temp_norm,2), size(temp_norm,3)]);
    [flimcoeff, scores, ~] = pca(vector_conv);
    filt_vector_PCA = (scores(:,1:filt_no_PC) * flimcoeff(:,1:filt_no_PC)') + mean(vector_conv, 1);
    filt_img_PCA = reshape(filt_vector_PCA, [size(temp_norm,1), size(temp_norm,2), size(temp_norm,3)]);
    stack1_norm = filt_img_PCA .* NoiseCorrection;
    stack1_norm = stack1_norm ./ max(stack1_norm(:));
    imwritestack(uint16(stack1_norm .* 65535), fullfile(out_dir, 'PCA_denoised.tif'));
else
    imwritestack(uint16(stack1_norm .* 65535), fullfile(out_dir, 'noPCA_input.tif'));
end

%% ====== Richardson-Lucy 光谱解混 ======
[H, W, ~] = size(stack1_norm);
N = H * W;

% --- 混合矩阵 M (F×2), 列 = 两条参考衰减曲线 ---
M = [ch1_g(1:frame), ch2_g(1:frame)];   % F×2
M = max(M, 0);                          % RL 要求非负
colsum = sum(M, 1)';                    % 2×1  ( = M' * ones(F,1), RL 归一化分母)

% --- 观测数据 D (F×N), 每列是一个像素的衰减曲线 ---
D = reshape(stack1_norm(:,:,1:frame), N, frame)';   % F×N
D = max(D, 0);

% --- 分块 RL 迭代 ---
C           = zeros(2, N);              % 2×N 丰度估计 (锁定在 num_iters 步)
keep_iters  = sort(save_iters);        % 捕获全部快照, 不受 num_iters 限制
K           = numel(keep_iters);
total_iters = max([num_iters, keep_iters]);   % 跑到足够远以捕获所有快照
iter_c1_flat = zeros(N, K);            % 各 save_iters 步的 c1 占比 (N×K)

fprintf('  RLSU: final=%d iters (跑到 %d 捕获全部快照), chunk=%d px ...\n', ...
    num_iters, total_iters, chunk_px);
tStart = tic;
for p0 = 1:chunk_px:N
    p1   = min(p0 + chunk_px - 1, N);
    Dblk = D(:, p0:p1);                 % F×n
    Cblk = ones(2, p1 - p0 + 1);        % 2×n
    ks = 1;
    for it = 1:total_iters
        Hu    = M * Cblk;               % F×n  正向模型
        ratio = Dblk ./ (Hu + eps);     % F×n  实测/模拟
        Cblk  = Cblk .* (M' * ratio) ./ colsum;   % 2×n  乘法更新 (colsum 广播)
        if it == num_iters                         % 锁定最终结果
            C(:, p0:p1) = Cblk;
        end
        if ks <= K && it == keep_iters(ks)         % 捕获收敛快照
            tblk = Cblk(1,:) + Cblk(2,:) + eps;
            iter_c1_flat(p0:p1, ks) = (Cblk(1,:) ./ tblk)';
            ks = ks + 1;
        end
    end
end
fprintf('  RLSU done in %.2f s\n', toc(tStart));

% --- 比例归一化 (sum-to-one), 与 fmincon/phasor 接口兼容 ---
c1 = C(1, :);  c2 = C(2, :);
tot = c1 + c2 + eps;
c1_map = (c1 ./ tot)';                  % N×1
c2_map = (c2 ./ tot)';

% --- Ratio-map TV 正则: 对 c2 做 TV, c1 由 1-c2 导出 (保 sum-to-one) ---
if tv_on
    c2_tv  = tv_denoise_rof(reshape(c2_map, H, W), tv_weight, tv_iter);
    c2_map = max(0, min(1, c2_tv(:)));
    c1_map = 1 - c2_map;
    fprintf('  [RLSU] ratio-map TV: weight=%.3f, iter=%d\n', tv_weight, tv_iter);
end

% --- 背景掩膜 (首帧强度低于阈值置零) ---
bg_mask = stack1_norm(:,:,1) < bg_thresh;     % H×W
c1_img  = reshape(c1_map, H, W);
c2_img  = reshape(c2_map, H, W);
c1_img(bg_mask) = 0;
c2_img(bg_mask) = 0;

% 组装 comp_raw (H×W×2)
comp_raw = cat(3, c1_img, c2_img);
amp_raw  = comp_raw;   % 与其他版本接口兼容

% --- 同时保存未归一化的原始丰度 (RLSU 的独有信息: 不强制 sum=1) ---
abund1 = reshape(c1', H, W);
abund2 = reshape(c2', H, W);

%% ====== 收敛序列 (每帧 = 一个 save_iters 步的 c1 占比) ======
cp   = 1 + paddingfactor;
crop = @(x) x(cp:end-paddingfactor, cp:end-paddingfactor, :);

if K > 0
    I1       = stack1_for_TMI(:,:,1);
    iter_raw = reshape(iter_c1_flat, H, W, K);   % 各步 c1 占比 (未 TV)

    % 待写出的序列集合: 第1组=原始 RL; tv_on 时追加第2组=TV 后, 文件名加 _TV
    src = {iter_raw, ''};
    if tv_on
        iter_tv = zeros(H, W, K);
        for k = 1:K
            ck_c2 = tv_denoise_rof(1 - iter_raw(:,:,k), tv_weight, tv_iter);
            iter_tv(:,:,k) = max(0, min(1, 1 - ck_c2));   % c1 = 1 - TV(c2)
        end
        src(2,:) = {iter_tv, '_TV'};
    end

    % 每组都写出 ratio_ch1 / TMIch1 / TMIch2 三个序列 (背景置零, 比例×首帧后逐帧归一化)
    for s = 1:size(src, 1)
        c1s = src{s,1};  suffix = src{s,2};
        c2s = 1 - c1s;
        T1  = zeros(H, W, K);  T2 = zeros(H, W, K);
        for k = 1:K
            a = c1s(:,:,k);  a(bg_mask) = 0;  c1s(:,:,k) = a;
            b = c2s(:,:,k);  b(bg_mask) = 0;  c2s(:,:,k) = b;
            u = a .* I1;  T1(:,:,k) = u ./ max(u(:) + eps);
            v = b .* I1;  T2(:,:,k) = v ./ max(v(:) + eps);
        end
        imwritestack(uint16(crop(c1s) .* 65535), fullfile(out_dir, ['ratio_ch1_iterations' suffix '.tif']));
        imwritestack(uint16(crop(T1)  .* 65535), fullfile(out_dir, ['TMIch1_iterations'    suffix '.tif']));
        imwritestack(uint16(crop(T2)  .* 65535), fullfile(out_dir, ['TMIch2_iterations'    suffix '.tif']));
    end
    fprintf('    收敛序列 (iters=%s) → ratio_ch1/TMIch1/TMIch2 _iterations(.|_TV).tif\n', mat2str(keep_iters));
end

%% ====== TMI 输出通道 ======
TMIch1 = comp_raw(:,:,1) .* stack1_for_TMI(:,:,1);
TMIch2 = comp_raw(:,:,2) .* stack1_for_TMI(:,:,1);
TMIch1 = TMIch1 ./ max(TMIch1(:) + eps);
TMIch2 = TMIch2 ./ max(TMIch2(:) + eps);

imwritestack(uint16(crop(TMIch1) .* 65535), fullfile(out_dir, 'TMIch1.tif'));
imwritestack(uint16(crop(TMIch2) .* 65535), fullfile(out_dir, 'TMIch2.tif'));
imwritestack(uint16(crop(comp_raw(:,:,1)) .* 65535), fullfile(out_dir, 'ratio_ch1.tif'));
imwritestack(uint16(crop(comp_raw(:,:,2)) .* 65535), fullfile(out_dir, 'ratio_ch2.tif'));
imwritestack(uint16(crop(abund1 ./ max(abund1(:)+eps)) .* 65535), fullfile(out_dir, 'abund_ch1.tif'));
imwritestack(uint16(crop(abund2 ./ max(abund2(:)+eps)) .* 65535), fullfile(out_dir, 'abund_ch2.tif'));

fprintf('  [step2_rlsu] 拆分完成 → %s\n', out_dir);
fprintf('    ratio ch1: [%.4f, %.4f]  ch2: [%.4f, %.4f]\n', ...
    min(c1_img,[],'all'), max(c1_img,[],'all'), ...
    min(c2_img,[],'all'), max(c2_img,[],'all'));

% =========================================================================
% 局部函数: ROF 全变分去噪 (Chambolle 2004 对偶投影算法)
% =========================================================================
function u = tv_denoise_rof(f, weight, n_iter)
% 求解 ROF 模型:   min_u  1/2 ||u - f||^2 + weight * TV(u)
%   TV(u) = sum sqrt(ux^2 + uy^2)  (各向同性全变分)
%   - 数据项拉住观测 f, TV 项压平内部噪声同时保锐边界
%   - weight (=λ) 越大越平滑
%   - Chambolle 对偶变量 p 投影迭代, tau<=1/8 保证收敛
[M, N] = size(f);
px = zeros(M, N);  py = zeros(M, N);
tau = 0.125;                       % <= 1/8
for k = 1:n_iter
    divp     = tv_div(px, py);
    [gx, gy] = tv_grad(divp - f / weight);
    gn = sqrt(gx.^2 + gy.^2);
    px = (px + tau * gx) ./ (1 + tau * gn);
    py = (py + tau * gy) ./ (1 + tau * gn);
end
u = f - weight * tv_div(px, py);
end

% --- 前向差分梯度 (Neumann 边界: 边缘处差分为 0) ---
function [gx, gy] = tv_grad(u)
[M, N] = size(u);
gx = zeros(M, N);  gy = zeros(M, N);
gx(:, 1:N-1) = u(:, 2:N) - u(:, 1:N-1);   % x 方向 (列)
gy(1:M-1, :) = u(2:M, :) - u(1:M-1, :);   % y 方向 (行)
end

% --- 散度 (前向梯度的负伴随, 满足 <grad u, p> = -<u, div p>) ---
function d = tv_div(px, py)
[M, N] = size(px);
dx = px;
dx(:, 2:N-1) = px(:, 2:N-1) - px(:, 1:N-2);
dx(:, N)     = -px(:, N-1);
dy = py;
dy(2:M-1, :) = py(2:M-1, :) - py(1:M-2, :);
dy(M, :)     = -py(M-1, :);
d = dx + dy;
end

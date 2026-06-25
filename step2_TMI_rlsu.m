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
pca_on     = 0;          % PCA 去噪开关 (0=关闭, 用于对比)
num_iters  = 100;        % RL 迭代次数 (2 组分收敛快, 50~200 即可)
chunk_px   = 50000;      % 每块处理的像素数 (控制内存; 越大越快越占内存)
save_iters = [1 5 20 50 100 500 1000];   % 额外保存这些迭代步的比例图, 便于观察收敛

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
C = zeros(2, N);                        % 2×N 丰度估计
keep_iters = save_iters(save_iters <= num_iters);
K = numel(keep_iters);
iter_c1_flat = zeros(N, K);             % 各 save_iters 步的 c1 占比 (N×K)

fprintf('  RLSU: %d iters, chunk=%d px ...\n', num_iters, chunk_px);
tStart = tic;
for p0 = 1:chunk_px:N
    p1   = min(p0 + chunk_px - 1, N);
    Dblk = D(:, p0:p1);                 % F×n
    Cblk = ones(2, p1 - p0 + 1);        % 2×n
    ks = 1;
    for it = 1:num_iters
        Hu    = M * Cblk;               % F×n  正向模型
        ratio = Dblk ./ (Hu + eps);     % F×n  实测/模拟
        Cblk  = Cblk .* (M' * ratio) ./ colsum;   % 2×n  乘法更新 (colsum 广播)
        if ks <= K && it == keep_iters(ks)         % 捕获收敛快照
            tblk = Cblk(1,:) + Cblk(2,:) + eps;
            iter_c1_flat(p0:p1, ks) = (Cblk(1,:) ./ tblk)';
            ks = ks + 1;
        end
    end
    C(:, p0:p1) = Cblk;
end
fprintf('  RLSU done in %.2f s\n', toc(tStart));

% --- 比例归一化 (sum-to-one), 与 fmincon/phasor 接口兼容 ---
c1 = C(1, :);  c2 = C(2, :);
tot = c1 + c2 + eps;
c1_map = (c1 ./ tot)';                  % N×1
c2_map = (c2 ./ tot)';

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
    iter_c1 = reshape(iter_c1_flat, H, W, K);
    for k = 1:K
        tmp = iter_c1(:,:,k);  tmp(bg_mask) = 0;  iter_c1(:,:,k) = tmp;
    end
    imwritestack(uint16(crop(iter_c1) .* 65535), ...
        fullfile(out_dir, 'ratio_ch1_iterations.tif'));
    fprintf('    收敛序列 (iters=%s) → ratio_ch1_iterations.tif\n', mat2str(keep_iters));
end

%% ====== TMI 输出通道 ======
TMIch1 = comp_raw(:,:,1) .* imgstack(:,:,1);
TMIch2 = comp_raw(:,:,2) .* imgstack(:,:,1);
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

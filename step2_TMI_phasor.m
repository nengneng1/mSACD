% =========================================================================
% step2_TMI_phasor.m — PCA去噪 + TMI 通道拆分 (Phasor 相量解混)
% =========================================================================
% 输入 (workspace): stack1_for_TMI, imgstack, frame, paddingfactor,
%                   FILE_OUT_DIR, fname, confname
% 输出 (workspace): amp_raw, comp_raw, TMIch1, TMIch2, H, W
% 输出 (tif):       {FILE_OUT_DIR}/TMI_phasor/
%
% 解混原理:
%   对每个像素的衰减曲线取第一谐波 Fourier 系数（相量坐标 g,s），
%   两条参考曲线 ch1/ch2 各对应 phasor 图上一个点，
%   混合像素落在两点连线上，投影位置即为混合比例 c1。
%   全程向量化，无迭代优化。
% =========================================================================

%% ====== 路径 ======
STEP_NAME = 'TMI_phasor';
addpath(genpath('F'));

%% ====== 参数 ======
gamma_on      = 1;
gamma_1       = 1;
bg_thresh     = 50/65535;
harmonics     = 1;      % 使用的谐波阶数（1=基频；可设为 [1,2] 使用多谐波投票）
phasor_filt   = 'median5';   % phasor 空间滤波: 'none'|'median3'|'median5'|'wiener3'

%% ====== 输出目录 ======
out_dir = fullfile(FILE_OUT_DIR, STEP_NAME);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ====== 衰减曲线 ======
load([confname, '.mat'], 'ch1', 'ch2');
ch1 = real(ch1(:));  ch2 = real(ch2(:));   % 确保列向量、实数

%% ====== Gamma 校正 ======
if gamma_on == 1
    stack1_g = stack1_for_TMI .^ gamma_1;
    ch1_g = (ch1 .^ gamma_1) ./ max(ch1 .^ gamma_1);
    ch2_g = (ch2 .^ gamma_1) ./ max(ch2 .^ gamma_1);
else
    stack1_g = stack1_for_TMI; ch1_g = ch1; ch2_g = ch2;
end

%% ====== 参考曲线可分性诊断 ======
R_ref = [ch1_g(1:frame), ch2_g(1:frame)];
c_ref = corrcoef(R_ref);
fprintf('  Reference: corr(ch1,ch2)=%.4f, cond([ch1,ch2])=%.2f\n', ...
    c_ref(1,2), cond(R_ref));

%% ====== PCA 去噪 ======
stack1_norm = double(stack1_g ./ max(stack1_g(:)));
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

%% ====== Phasor 相量解混 ======
[H, W, ~] = size(stack1_norm);
t = (0:frame-1)';   % 时间轴（帧序号）

% --- 参考曲线 phasor 坐标（每个谐波一个点） ---
% 存为 [n_harmonics × 2]，列: [g, s]
ref1 = zeros(numel(harmonics), 2);
ref2 = zeros(numel(harmonics), 2);
for hi = 1:numel(harmonics)
    n = harmonics(hi);
    omega = 2*pi*n / frame;
    c1n = ch1_g(1:frame);  c2n = ch2_g(1:frame);
    ref1(hi,:) = [dot(c1n, cos(omega*t)), dot(c1n, sin(omega*t))] / sum(c1n);
    ref2(hi,:) = [dot(c2n, cos(omega*t)), dot(c2n, sin(omega*t))] / sum(c2n);
end

fprintf('  Phasor ref1: ');
fprintf('n=%d (g=%.4f, s=%.4f)  ', [harmonics; ref1(:,1)'; ref1(:,2)']);
fprintf('\n  Phasor ref2: ');
fprintf('n=%d (g=%.4f, s=%.4f)  ', [harmonics; ref2(:,1)'; ref2(:,2)']);
fprintf('\n');

% --- 每像素 phasor 坐标（向量化）---
% stack1_norm: H×W×frame → 展开为 (H*W)×frame
S = reshape(stack1_norm(:,:,1:frame), H*W, frame);   % (HW)×F
I_sum = sum(S, 2) + eps;                              % (HW)×1

% 对每个谐波计算 g, s，取平均投影
alpha_accum = zeros(H*W, 1);
omega1 = 2*pi*harmonics(1) / frame;
g_raw  = (S * cos(omega1*t)) ./ I_sum;   % 保留未滤波版本用于 phasor 图对比
s_raw  = (S * sin(omega1*t)) ./ I_sum;

for hi = 1:numel(harmonics)
    n = harmonics(hi);
    omega = 2*pi*n / frame;
    g_px = (S * cos(omega*t)) ./ I_sum;   % (HW)×1
    s_px = (S * sin(omega*t)) ./ I_sum;   % (HW)×1

    % --- Phasor 空间滤波 ---
    g_img = reshape(g_px, H, W);
    s_img = reshape(s_px, H, W);
    switch phasor_filt
        case 'median3'
            g_img = medfilt2(g_img, [3 3]);
            s_img = medfilt2(s_img, [3 3]);
        case 'median5'
            g_img = medfilt2(g_img, [5 5]);
            s_img = medfilt2(s_img, [5 5]);
        case 'wiener3'
            g_img = wiener2(g_img, [3 3]);
            s_img = wiener2(s_img, [3 3]);
        % case 'none': 不做任何处理
    end
    g_px = g_img(:);
    s_px = s_img(:);

    % 投影到 ref1–ref2 连线
    dg = ref1(hi,1) - ref2(hi,1);
    ds = ref1(hi,2) - ref2(hi,2);
    denom = dg^2 + ds^2 + eps;
    c1_map = ((g_px - ref2(hi,1))*dg + (s_px - ref2(hi,2))*ds) / denom;
    alpha_accum = alpha_accum + c1_map;
end
g_filt = g_px;   % 保留已滤波版本（harmonics(1) 对应的最后一次循环值）
s_filt = s_px;
c1_map = alpha_accum / numel(harmonics);   % 多谐波平均
c1_map = max(0, min(1, c1_map));           % 截断到 [0,1]
c2_map = 1 - c1_map;

% --- 背景掩膜（首帧强度低于阈值的像素置零）---
bg_mask = stack1_norm(:,:,1) < bg_thresh;   % H×W
c1_img  = reshape(c1_map, H, W);
c2_img  = reshape(c2_map, H, W);
c1_img(bg_mask) = 0;
c2_img(bg_mask) = 0;

% 组装 comp_raw (H×W×2)
comp_raw = cat(3, c1_img, c2_img);
amp_raw  = comp_raw;   % 与 fmincon 版本接口兼容

%% ====== 保存 phasor 图（滤波前 + 滤波后对比）======
n_samp    = min(5000, H*W);
idx_s     = randperm(H*W, n_samp);
th        = linspace(0, pi, 200);
phasor_ax = [-0.1 1.1 -0.1 0.7];   % [xmin xmax ymin ymax]

% --- 滤波前 ---
fig = figure('Visible','off');
scatter(g_raw(idx_s), s_raw(idx_s), 1, [0.7 0.7 0.7]); hold on;
plot(0.5+0.5*cos(th), 0.5*sin(th), 'k--', 'LineWidth', 0.8);
plot(ref1(1,1), ref1(1,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName','ch1');
plot(ref2(1,1), ref2(1,2), 'bo', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName','ch2');
line([ref1(1,1) ref2(1,1)], [ref1(1,2) ref2(1,2)], 'Color','k');
legend; xlabel('g'); ylabel('s'); axis equal; axis(phasor_ax);
title(sprintf('Phasor raw  n=%d  |  %s', harmonics(1), fname), 'Interpreter','none');
saveas(fig, fullfile(out_dir, 'phasor_raw.png'));
close(fig);

% --- 滤波后 ---
fig = figure('Visible','off');
scatter(g_filt(idx_s), s_filt(idx_s), 1, [0.7 0.7 0.7]); hold on;
plot(0.5+0.5*cos(th), 0.5*sin(th), 'k--', 'LineWidth', 0.8);
plot(ref1(1,1), ref1(1,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName','ch1');
plot(ref2(1,1), ref2(1,2), 'bo', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName','ch2');
line([ref1(1,1) ref2(1,1)], [ref1(1,2) ref2(1,2)], 'Color','k');
legend; xlabel('g'); ylabel('s'); axis equal; axis(phasor_ax);
title(sprintf('Phasor %s  n=%d  |  %s', phasor_filt, harmonics(1), fname), 'Interpreter','none');
saveas(fig, fullfile(out_dir, 'phasor_filt.png'));
close(fig);

%% ====== TMI 输出通道 ======
TMIch1 = comp_raw(:,:,1) .* imgstack(:,:,1);
TMIch2 = comp_raw(:,:,2) .* imgstack(:,:,1);
TMIch1 = TMIch1 ./ max(TMIch1(:) + eps);
TMIch2 = TMIch2 ./ max(TMIch2(:) + eps);

cp   = 1 + paddingfactor;
crop = @(x) x(cp:end-paddingfactor, cp:end-paddingfactor, :);

imwritestack(uint16(crop(TMIch1) .* 65535), fullfile(out_dir, 'TMIch1.tif'));
imwritestack(uint16(crop(TMIch2) .* 65535), fullfile(out_dir, 'TMIch2.tif'));
imwritestack(uint16(crop(comp_raw(:,:,1)) .* 65535), fullfile(out_dir, 'ratio_ch1.tif'));
imwritestack(uint16(crop(comp_raw(:,:,2)) .* 65535), fullfile(out_dir, 'ratio_ch2.tif'));

fprintf('  [step2_phasor] 拆分完成 → %s\n', out_dir);
fprintf('    comp ch1: [%.4f, %.4f]  ch2: [%.4f, %.4f]\n', ...
    min(c1_img,[],'all'), max(c1_img,[],'all'), ...
    min(c2_img,[],'all'), max(c2_img,[],'all'));

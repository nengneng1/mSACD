% =========================================================================
% step2_TMI_phasor.m — PCA去噪 + TMI 通道拆分 (Phasor 相量解混)
% =========================================================================
% 输入 (workspace): stack1_for_TMI, imgstack, frame, paddingfactor,
%                   FILE_OUT_DIR, fname, confname
% 输出 (workspace): comp_raw, TMIch1, TMIch2, H, W
% 输出 (tif/png):   {FILE_OUT_DIR}/TMI_phasor/
%
% 解混原理:
%   对每像素的衰减曲线取第一谐波 Fourier 系数 (g,s)，
%   两条参考曲线 ch1/ch2 各映射为 phasor 图上一个点 (ref1, ref2)，
%   混合像素落在 ref2→ref1 连线上，投影位置即为 c1（ch1 占比）。
%   全程向量化，无迭代优化。
% =========================================================================

%% ====== 路径与标签 ======
STEP_NAME = 'TMI_phasor';
UNMIX_TAG = 'phasor';    % 供 step2b / step3 区分输出子目录
addpath(genpath('F'));

%% ====== 参数 ======
% --- 通用 ---
bg_thresh = 0/65535;
pca_on    = 1;          % PCA 去噪: 1=开启, 0=关闭

% --- Phasor ---
harmonics   = [1];      % 谐波阶数; [1] 基频; [1,2] 多谐波均值投影
phasor_filt = 'cwf';    % 相量空间滤波: 'none'|'median3'|'median5'|'wiener3'|'cwf'

%% ====== 输出目录 ======
out_dir = fullfile(FILE_OUT_DIR, STEP_NAME);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ====== 加载参考曲线 ======
load([confname, '.mat'], 'ch1', 'ch2');
ch1 = real(ch1(:));
ch2 = real(ch2(:));

%% ====== 参考曲线 (无 gamma, 保持线性混合假设) ======
stack1_g = stack1_for_TMI;
ch1_g    = ch1;
ch2_g    = ch2;

%% ====== 参考曲线可分性诊断 ======
R_ref = [ch1_g(1:frame), ch2_g(1:frame)];
c_ref = corrcoef(R_ref);
fprintf('  Reference: corr(ch1,ch2)=%.4f, cond([ch1,ch2])=%.2f\n', ...
    c_ref(1,2), cond(R_ref));

%% ====== PCA 去噪 ======
stack1_norm = double(stack1_g ./ max(stack1_g(:)));
if pca_on
    filt_no_PC       = 3;
    NoiseCorrection  = sqrt(mean(stack1_norm, [1,2]) + eps);
    temp_norm        = stack1_norm ./ NoiseCorrection;
    temp_norm(temp_norm < 0) = 0;
    vector_conv      = reshape(temp_norm, [], size(temp_norm,3));
    [flimcoeff, scores, ~] = pca(vector_conv);
    filt_vector_PCA  = (scores(:,1:filt_no_PC) * flimcoeff(:,1:filt_no_PC)') ...
                       + mean(vector_conv, 1);
    filt_img_PCA     = reshape(filt_vector_PCA, size(temp_norm));
    stack1_norm      = filt_img_PCA .* NoiseCorrection;
    stack1_norm      = stack1_norm ./ max(stack1_norm(:));
    imwritestack(uint16(stack1_norm .* 65535), fullfile(out_dir, 'PCA_denoised.tif'));
else
    imwritestack(uint16(stack1_norm .* 65535), fullfile(out_dir, 'noPCA_input.tif'));
end

% 注: 解混前不再做 Dark sectioning。背景去除统一放到 step3 时变分离之后。

%% ====== Phasor：参考曲线坐标 ======
[H, W, ~] = size(stack1_norm);
t = (0:frame-1)';

ref1 = zeros(numel(harmonics), 2);
ref2 = zeros(numel(harmonics), 2);
for hi = 1:numel(harmonics)
    n     = harmonics(hi);
    omega = 2*pi*n / frame;
    c1n   = ch1_g(1:frame);
    c2n   = ch2_g(1:frame);
    d1    = sum(c1n);  d2 = sum(c2n);
    if d1 > eps, ref1(hi,:) = [dot(c1n, cos(omega*t)), dot(c1n, sin(omega*t))] / d1; end
    if d2 > eps, ref2(hi,:) = [dot(c2n, cos(omega*t)), dot(c2n, sin(omega*t))] / d2; end
end
fprintf('  Phasor ref1: ');
fprintf('n=%d (g=%.4f, s=%.4f)  ', [harmonics; ref1(:,1)'; ref1(:,2)']);
fprintf('\n  Phasor ref2: ');
fprintf('n=%d (g=%.4f, s=%.4f)  ', [harmonics; ref2(:,1)'; ref2(:,2)']);
fprintf('\n');

%% ====== Phasor：逐像素 g/s + 空间滤波 + 投影 ======
S      = reshape(stack1_norm(:,:,1:frame), H*W, frame);   % (HW)×F
I_sum  = sum(S, 2) + eps;

% 保留未滤波版本（仅第一谐波，供诊断图对比）
omega1 = 2*pi*harmonics(1) / frame;
g_raw  = (S * cos(omega1*t)) ./ I_sum;
s_raw  = (S * sin(omega1*t)) ./ I_sum;

alpha_accum = zeros(H*W, 1);

for hi = 1:numel(harmonics)
    n     = harmonics(hi);
    omega = 2*pi*n / frame;
    g_px  = (S * cos(omega*t)) ./ I_sum;
    s_px  = (S * sin(omega*t)) ./ I_sum;

    % --- 相量空间滤波 ---
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
        case 'cwf'
            % CWF：在未归一化 Fourier 系数上做 Anscombe + 小波去噪，
            % 之后再归一化，避免低强度像素中 I_sum 噪声放大 g/s 误差
            I_vec   = I_sum - eps;
            Fr_img  = reshape(g_px .* I_vec, H, W);
            Fi_img  = reshape(s_px .* I_vec, H, W);
            I_img   = reshape(I_vec,          H, W);
            ansc    = @(x) 2*sqrt(max(x,0) + 3/8);
            iansc   = @(y) max((y/2).^2 - 3/8, 0);
            if exist('wdenoise2', 'file')
                Fr_d = wdenoise2(ansc(Fr_img), 'Wavelet', 'db4');
                Fi_d = wdenoise2(ansc(Fi_img), 'Wavelet', 'db4');
                I_d  = wdenoise2(ansc(I_img),  'Wavelet', 'db4');
            else
                warning('[cwf] wdenoise2 未找到，回退到 Anscombe+medfilt2');
                Fr_d = medfilt2(ansc(Fr_img), [3 3]);
                Fi_d = medfilt2(ansc(Fi_img), [3 3]);
                I_d  = medfilt2(ansc(I_img),  [3 3]);
            end
            I_f   = max(iansc(I_d), eps);
            g_img = iansc(Fr_d) ./ I_f;
            s_img = iansc(Fi_d) ./ I_f;
            % 截断极端值（背景像素 I_f≈eps 会产生极大 g/s）
            g_img = max(min(g_img,  1.5), -0.3);
            s_img = max(min(s_img,  1.0), -0.3);
        % case 'none': 不处理
    end

    g_px = g_img(:);
    s_px = s_img(:);

    % 投影到 ref2→ref1 连线，得到 c1
    dg    = ref1(hi,1) - ref2(hi,1);
    ds    = ref1(hi,2) - ref2(hi,2);
    denom = dg^2 + ds^2 + eps;
    c1_map = ((g_px - ref2(hi,1))*dg + (s_px - ref2(hi,2))*ds) / denom;
    alpha_accum = alpha_accum + c1_map;
end

g_filt = g_px;   % 保留滤波后版本（最后一次循环值，对应 harmonics(1)）
s_filt = s_px;

c1_map = alpha_accum / numel(harmonics);
c1_map = max(0, min(1, c1_map));
c2_map = 1 - c1_map;

%% ====== 背景掩膜 + comp_raw 组装 ======
bg_mask = stack1_norm(:,:,1) < bg_thresh;
c1_img  = reshape(c1_map, H, W);
c2_img  = reshape(c2_map, H, W);
c1_img(bg_mask) = 0;
c2_img(bg_mask) = 0;

comp_raw = cat(3, c1_img, c2_img);

%% ====== Phasor 诊断图（重写版：简洁，散点必可见）======
% 前景像素 g/s：raw=未滤波, filt=空间滤波后；用 c1 给滤波散点着色
fg_idx = ~reshape(bg_mask, H*W, 1);
G_r = g_raw(fg_idx);    S_r = s_raw(fg_idx);
G_f = g_filt(fg_idx);   S_f = s_filt(fg_idx);
c1_fg = c1_map(fg_idx);
mu_c1 = mean(c1_fg);  sd_c1 = std(c1_fg);

% 随机下采样, 避免散点过密
rng(42);
n_pts  = numel(G_r);
n_samp = min(15000, n_pts);
samp   = randperm(n_pts, n_samp);

% 参考点 / 标定线 / 通用半圆
P1 = ref1(1,:);  P2 = ref2(1,:);            % P1=ch1, P2=ch2
tt   = linspace(-0.15, 1.15, 200);
calg = P2(1) + tt*(P1(1)-P2(1));  cals = P2(2) + tt*(P1(2)-P2(2));
ph   = linspace(0, pi, 200);
semg = 0.5 + 0.5*cos(ph);  sems = 0.5*sin(ph);

% 坐标范围：直接由"实际要画的散点 + 两个参考点"取 min/max + 8% 边距
% → 不再用任何写死的框, 散点一定落在画面内
allg = [G_r(samp); G_f(samp); P1(1); P2(1)];
alls = [S_r(samp); S_f(samp); P1(2); P2(2)];
mg = 0.08 * max(0.05, max(allg)-min(allg));
ms = 0.08 * max(0.05, max(alls)-min(alls));
ax = [min(allg)-mg, max(allg)+mg, min(alls)-ms, max(alls)+ms];

fprintf('  [Phasor] ref1=(%.3f,%.3f)  ref2=(%.3f,%.3f)  c1: mu=%.3f sd=%.3f\n', ...
    P1(1),P1(2), P2(1),P2(2), mu_c1, sd_c1);
fprintf('  [Phasor] 散点范围 g=[%.3f,%.3f]  s=[%.3f,%.3f]\n', ...
    min(allg),max(allg), min(alls),max(alls));

% ---- 图1: g/s 散点 (左=未滤波  右=滤波后, 按 c1 着色) ----
fig = figure('Visible','off', 'Position',[60 60 1300 560]);

subplot(1,2,1);
scatter(G_r(samp), S_r(samp), 14, [0.40 0.40 0.40], 'filled', 'MarkerFaceAlpha',0.40, 'HandleVisibility','off');
hold on;
plot(semg, sems, ':',  'Color',[0.6 0.6 0.6], 'LineWidth',1, 'HandleVisibility','off');
plot(calg, cals, 'k--', 'LineWidth',1.5, 'DisplayName','标定线 ref2→ref1');
plot(P1(1),P1(2), 'o', 'MarkerFaceColor',[0.90 0.20 0.20], 'MarkerEdgeColor','k', 'MarkerSize',12, 'DisplayName','ref1 (ch1)');
plot(P2(1),P2(2), 's', 'MarkerFaceColor',[0.20 0.50 0.95], 'MarkerEdgeColor','k', 'MarkerSize',12, 'DisplayName','ref2 (ch2)');
axis equal; axis(ax); grid on; box on;
xlabel('g'); ylabel('s'); legend('Location','best','FontSize',8);
title(sprintf('未滤波  (N=%d)', n_samp), 'FontSize',10);

subplot(1,2,2);
scatter(G_f(samp), S_f(samp), 14, c1_fg(samp), 'filled', 'MarkerFaceAlpha',0.60, 'HandleVisibility','off');
colormap(parula); caxis([0 1]); cb = colorbar; cb.Label.String = 'c_1 (ch1 占比)';
hold on;
plot(semg, sems, ':',  'Color',[0.6 0.6 0.6], 'LineWidth',1, 'HandleVisibility','off');
plot(calg, cals, 'k--', 'LineWidth',1.5, 'DisplayName','标定线 ref2→ref1');
plot(P1(1),P1(2), 'o', 'MarkerFaceColor',[0.90 0.20 0.20], 'MarkerEdgeColor','k', 'MarkerSize',12, 'DisplayName','ref1 (ch1)');
plot(P2(1),P2(2), 's', 'MarkerFaceColor',[0.20 0.50 0.95], 'MarkerEdgeColor','k', 'MarkerSize',12, 'DisplayName','ref2 (ch2)');
axis equal; axis(ax); grid on; box on;
xlabel('g'); ylabel('s'); legend('Location','best','FontSize',8);
title(sprintf('%s 滤波', phasor_filt), 'FontSize',10);

sgtitle(sprintf('Phasor  n=%d  |  %s', harmonics(1), fname), 'Interpreter','none', 'FontSize',11);
saveas(fig, fullfile(out_dir, 'phasor_scatter.png')); close(fig);

% ---- 图2: c1 分布直方图 ----
fig = figure('Visible','off', 'Position',[60 60 640 460]);
histogram(c1_fg, 60, 'FaceColor',[0.18 0.55 0.85], 'EdgeColor','none', 'Normalization','probability');
hold on;
xline(mu_c1, 'r-',  'LineWidth',2,   'Label',sprintf('μ=%.3f',mu_c1), 'FontSize',9);
xline(0.5,   'k--', 'LineWidth',1.2, 'Label','0.5', 'FontSize',9);
xlabel('c_1 (ch1 占比)', 'FontSize',11); ylabel('概率', 'FontSize',11);
xlim([0 1]); grid on; box on;
title(sprintf('解混比例分布  μ=%.3f  σ=%.3f  |  %s', mu_c1, sd_c1, fname), ...
    'Interpreter','none', 'FontSize',10);
saveas(fig, fullfile(out_dir, 'c1_distribution.png')); close(fig);

fprintf('  [Phasor 诊断图] → %s\n', out_dir);

%% ====== TMI 输出通道 ======
cp   = 1 + paddingfactor;
crop = @(x) x(cp:end-paddingfactor, cp:end-paddingfactor, :);

I_ref  = stack1_norm(:,:,1);      % 去背景+PCA后的第一帧作为强度权重
TMIch1 = comp_raw(:,:,1) .* I_ref;
TMIch1 = TMIch1 ./ max(TMIch1(:) + eps);
TMIch2 = comp_raw(:,:,2) .* I_ref;
TMIch2 = TMIch2 ./ max(TMIch2(:) + eps);

imwritestack(uint16(crop(TMIch1) .* 65535),           fullfile(out_dir, 'TMIch1.tif'));
imwritestack(uint16(crop(TMIch2) .* 65535),           fullfile(out_dir, 'TMIch2.tif'));
imwritestack(uint16(crop(comp_raw(:,:,1)) .* 65535),  fullfile(out_dir, 'ratio_ch1.tif'));
imwritestack(uint16(crop(comp_raw(:,:,2)) .* 65535),  fullfile(out_dir, 'ratio_ch2.tif'));

fprintf('  [step2_phasor] 拆分完成 → %s\n', out_dir);
fprintf('    comp ch1: [%.4f, %.4f]  ch2: [%.4f, %.4f]\n', ...
    min(c1_img,[],'all'), max(c1_img,[],'all'), ...
    min(c2_img,[],'all'), max(c2_img,[],'all'));

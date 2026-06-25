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
phasor_filt = 'median3';    % 相量空间滤波: 'none'|'median3'|'median5'|'wiener3'|'cwf'

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
    filt_no_PC       = 2;
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

%% ====== Phasor 诊断图 ======
rng(42);
fg_idx    = ~reshape(bg_mask, H*W, 1);
fg_all    = find(fg_idx);

G_r = g_raw(fg_idx);    S_r = s_raw(fg_idx);
G_f = g_filt(fg_idx);   S_f = s_filt(fg_idx);
gc_r = mean(G_r);  sc_r = mean(S_r);
gc_f = mean(G_f);  sc_f = mean(S_f);

% 自适应轴范围（数据 2%~98% 分位 + 参考点，外扩 25%）
g_ext = [G_r; G_f; ref1(1,1); ref2(1,1)];
s_ext = [S_r; S_f; ref1(1,2); ref2(1,2)];
gm    = 0.25 * max(0.05, prctile(g_ext,98) - prctile(g_ext,2));
sm    = 0.25 * max(0.05, prctile(s_ext,98) - prctile(s_ext,2));
phasor_zoom = [prctile(g_ext,2)-gm, prctile(g_ext,98)+gm, ...
               prctile(s_ext,2)-sm, prctile(s_ext,98)+sm];
phasor_full = [-0.1 1.1 -0.1 0.7];

% 数据分布主轴（SVD，不依赖 Statistics Toolbox）
[~,~,Vr] = svd([G_r(:)-gc_r, S_r(:)-sc_r], 'econ');  dir_r = Vr(:,1);
[~,~,Vf] = svd([G_f(:)-gc_f, S_f(:)-sc_f], 'econ');  dir_f = Vf(:,1);

% 主轴方向对齐 ref2→ref1（避免 180° 翻转）
d_ref  = [ref1(1,1)-ref2(1,1), ref1(1,2)-ref2(1,2)];
len_ref = norm(d_ref);
d_unit  = d_ref / len_ref;
if dot(dir_r', d_ref) < 0, dir_r = -dir_r; end
if dot(dir_f', d_ref) < 0, dir_f = -dir_f; end

% 拟合直线（主轴方向，±0.5 延伸）
t_ext   = linspace(-0.5, 0.5, 300);
fit_g_r = gc_r + t_ext*dir_r(1);  fit_s_r = sc_r + t_ext*dir_r(2);
fit_g_f = gc_f + t_ext*dir_f(1);  fit_s_f = sc_f + t_ext*dir_f(2);

% 标定线（ref2→ref1，两端各延伸 20%）
t_cal = linspace(-0.2, 1.2, 300);
cal_g = ref2(1,1) + t_cal*d_ref(1);
cal_s = ref2(1,2) + t_cal*d_ref(2);

% 半圆曲线（phasor 通用背景）
th = linspace(0, pi, 300);

% 定量偏差：质心到标定线的有符号垂直距离和主轴角偏差
nvec  = [-d_ref(2), d_ref(1)] / len_ref;
off_r = dot([gc_r-ref2(1,1), sc_r-ref2(1,2)], nvec);
off_f = dot([gc_f-ref2(1,1), sc_f-ref2(1,2)], nvec);

ang_cal = atan2d(d_ref(2), d_ref(1));
dang_r  = atan2d(dir_r(2), dir_r(1)) - ang_cal;
dang_f  = atan2d(dir_f(2), dir_f(1)) - ang_cal;
if abs(dang_r) > 90, dang_r = dang_r - sign(dang_r)*180; end
if abs(dang_f) > 90, dang_f = dang_f - sign(dang_f)*180; end

% 质心在标定线上的垂足（用于偏移箭头）
t_r    = dot([gc_r-ref2(1,1), sc_r-ref2(1,2)], d_unit);
proj_r = [ref2(1,1)+t_r*d_unit(1), ref2(1,2)+t_r*d_unit(2)];
t_f    = dot([gc_f-ref2(1,1), sc_f-ref2(1,2)], d_unit);
proj_f = [ref2(1,1)+t_f*d_unit(1), ref2(1,2)+t_f*d_unit(2)];

fprintf('  [Phasor 诊断]\n');
fprintf('    质心偏移(⊥标定线): raw=%+.4f  filt=%+.4f\n', off_r, off_f);
fprintf('    主轴角偏差:         raw=%+.2f°  filt=%+.2f°\n', dang_r, dang_f);
fprintf('    前景质心 raw  = (%.4f, %.4f)\n', gc_r, sc_r);
fprintf('    前景质心 filt = (%.4f, %.4f)\n', gc_f, sc_f);

% --- 散点采样（随机取子集，避免过密）---
n_fg_samp = min(8000, numel(fg_all));
samp_idx  = randperm(numel(fg_all), n_fg_samp);
G_rs = G_r(samp_idx);   S_rs = S_r(samp_idx);
G_fs = G_f(samp_idx);   S_fs = S_f(samp_idx);
c1_s = c1_map(fg_all(samp_idx));   % 按 c1 着色

% ---- 图1: RAW 散点图（全局）----
fig1 = figure('Visible','off', 'Position',[50 50 760 640]);
scatter(G_rs, S_rs, 4, [0.55 0.55 0.55], 'filled', 'MarkerFaceAlpha',0.25, 'HandleVisibility','off');
hold on;
plot(0.5+0.5*cos(th), 0.5*sin(th), 'Color',[0.7 0.7 0.7], 'LineStyle','-.', 'LineWidth',1.2, 'HandleVisibility','off');
plot(cal_g, cal_s, 'k--', 'LineWidth',2,   'DisplayName','标定线 (ref1-ref2)');
plot(fit_g_r, fit_s_r, 'Color',[0 0.75 0.75], 'LineWidth',2.5, ...
    'DisplayName', sprintf('数据主轴  offset=%+.4f  Δθ=%+.2f°', off_r, dang_r));
plot(ref1(1,1), ref1(1,2), 'o', 'Color','k', 'MarkerFaceColor',[1 0.3 0.3], 'MarkerSize',14, 'LineWidth',2, 'DisplayName','ref1 (ch1)');
plot(ref2(1,1), ref2(1,2), 's', 'Color','k', 'MarkerFaceColor',[0.3 0.6 1], 'MarkerSize',14, 'LineWidth',2, 'DisplayName','ref2 (ch2)');
plot(gc_r, sc_r, '^', 'Color','k', 'MarkerFaceColor',[0 0.75 0.75], 'MarkerSize',12, 'LineWidth',1.5, ...
    'DisplayName', sprintf('质心 (%.4f, %.4f)', gc_r, sc_r));
quiver(proj_r(1), proj_r(2), gc_r-proj_r(1), sc_r-proj_r(2), 0, ...
    'Color',[0 0.75 0.75], 'LineWidth',2.5, 'MaxHeadSize',0.8, 'HandleVisibility','off');
legend('Location','northeast', 'FontSize',8);
xlabel('g (cos 分量)', 'FontSize',11); ylabel('s (sin 分量)', 'FontSize',11);
axis equal; axis(phasor_full); box on;
title(sprintf('Phasor [RAW]  n=%d  |  %s', harmonics(1), fname), 'Interpreter','none', 'FontSize',10);
saveas(fig1, fullfile(out_dir, 'phasor_diag_raw.png')); close(fig1);

% ---- 图2: 滤波后散点图（按 c1 着色，全局）----
fig2 = figure('Visible','off', 'Position',[50 50 760 640]);
scatter(G_fs, S_fs, 4, c1_s, 'filled', 'MarkerFaceAlpha',0.35, 'HandleVisibility','off');
colormap(parula); hcb = colorbar; hcb.Label.String = 'c₁ (ch1 占比)';
hold on;
plot(0.5+0.5*cos(th), 0.5*sin(th), 'Color',[0.7 0.7 0.7], 'LineStyle','-.', 'LineWidth',1.2, 'HandleVisibility','off');
plot(cal_g, cal_s, 'k--', 'LineWidth',2,   'DisplayName','标定线 (ref1-ref2)');
plot(fit_g_f, fit_s_f, 'Color',[0.85 0.65 0], 'LineWidth',2.5, ...
    'DisplayName', sprintf('数据主轴  offset=%+.4f  Δθ=%+.2f°', off_f, dang_f));
plot(ref1(1,1), ref1(1,2), 'o', 'Color','k', 'MarkerFaceColor',[1 0.3 0.3], 'MarkerSize',14, 'LineWidth',2, 'DisplayName','ref1 (ch1)');
plot(ref2(1,1), ref2(1,2), 's', 'Color','k', 'MarkerFaceColor',[0.3 0.6 1], 'MarkerSize',14, 'LineWidth',2, 'DisplayName','ref2 (ch2)');
plot(gc_f, sc_f, '^', 'Color','k', 'MarkerFaceColor',[0.85 0.65 0], 'MarkerSize',12, 'LineWidth',1.5, ...
    'DisplayName', sprintf('质心 (%.4f, %.4f)', gc_f, sc_f));
quiver(proj_f(1), proj_f(2), gc_f-proj_f(1), sc_f-proj_f(2), 0, ...
    'Color',[0.85 0.65 0], 'LineWidth',2.5, 'MaxHeadSize',0.8, 'HandleVisibility','off');
legend('Location','northeast', 'FontSize',8);
xlabel('g (cos 分量)', 'FontSize',11); ylabel('s (sin 分量)', 'FontSize',11);
axis equal; axis(phasor_full); box on;
title(sprintf('Phasor [%s]  n=%d  |  %s', phasor_filt, harmonics(1), fname), 'Interpreter','none', 'FontSize',10);
saveas(fig2, fullfile(out_dir, 'phasor_diag_filt.png')); close(fig2);

% ---- 图3: RAW vs 滤波 对比（左=RAW缩放  右=滤波全局）----
fig3 = figure('Visible','off', 'Position',[50 50 1400 600]);

subplot(1,2,1);
scatter(G_rs, S_rs, 4, [0.55 0.55 0.55], 'filled', 'MarkerFaceAlpha',0.25, 'HandleVisibility','off');
hold on;
plot(cal_g, cal_s, 'k--', 'LineWidth',1.8);
plot(fit_g_r, fit_s_r, 'Color',[0 0.75 0.75], 'LineWidth',2.2);
plot(ref1(1,1), ref1(1,2), 'o', 'Color','k', 'MarkerFaceColor',[1 0.3 0.3], 'MarkerSize',12, 'LineWidth',2);
plot(ref2(1,1), ref2(1,2), 's', 'Color','k', 'MarkerFaceColor',[0.3 0.6 1], 'MarkerSize',12, 'LineWidth',2);
quiver(proj_r(1),proj_r(2),gc_r-proj_r(1),sc_r-proj_r(2),0,'Color',[0 0.75 0.75],'LineWidth',2,'MaxHeadSize',1);
axis equal; axis(phasor_zoom); box on;
xlabel('g'); ylabel('s');
title(sprintf('RAW [缩放]  offset=%+.4f  Δθ=%+.1f°', off_r, dang_r), 'FontSize',10);

subplot(1,2,2);
scatter(G_fs, S_fs, 4, c1_s, 'filled', 'MarkerFaceAlpha',0.35, 'HandleVisibility','off');
colormap(parula); colorbar;
hold on;
plot(0.5+0.5*cos(th), 0.5*sin(th), 'Color',[0.7 0.7 0.7], 'LineStyle','-.', 'LineWidth',1);
plot(cal_g, cal_s, 'k--', 'LineWidth',1.8, 'DisplayName','标定线');
plot(fit_g_f, fit_s_f, 'Color',[0.85 0.65 0], 'LineWidth',2.2, 'DisplayName','数据主轴');
plot(ref1(1,1), ref1(1,2), 'o', 'Color','k', 'MarkerFaceColor',[1 0.3 0.3], 'MarkerSize',12, 'LineWidth',2, 'DisplayName','ref1');
plot(ref2(1,1), ref2(1,2), 's', 'Color','k', 'MarkerFaceColor',[0.3 0.6 1], 'MarkerSize',12, 'LineWidth',2, 'DisplayName','ref2');
quiver(proj_f(1),proj_f(2),gc_f-proj_f(1),sc_f-proj_f(2),0,'Color',[0.85 0.65 0],'LineWidth',2,'MaxHeadSize',1);
legend('Location','northeast', 'FontSize',8);
axis equal; axis(phasor_full); box on;
xlabel('g'); ylabel('s');
title(sprintf('%s [全局]  offset=%+.4f  Δθ=%+.1f°', phasor_filt, off_f, dang_f), 'FontSize',10);

sgtitle(sprintf('Phasor 散点对比  n=%d  |  %s', harmonics(1), fname), 'Interpreter','none', 'FontSize',11);
saveas(fig3, fullfile(out_dir, 'phasor_compare.png')); close(fig3);

% ---- 图4: c1 比例分布直方图 ----
fig4 = figure('Visible','off', 'Position',[50 50 640 460]);
c1_fg = c1_map(fg_idx);
histogram(c1_fg, 60, 'FaceColor',[0.18 0.55 0.85], 'EdgeColor','none', ...
    'Normalization','probability');
hold on;
xline(mean(c1_fg), 'r-', 'LineWidth',2, ...
    'Label', sprintf('均值=%.3f', mean(c1_fg)), 'LabelVerticalAlignment','bottom', 'FontSize',9);
xline(0.5, 'k--', 'LineWidth',1.5, ...
    'Label', '0.5', 'LabelVerticalAlignment','bottom', 'FontSize',9);
xlabel('c₁（ch1 占比）', 'FontSize',11);
ylabel('概率', 'FontSize',11);
title(sprintf('解混比例分布  μ=%.3f  σ=%.3f  |  %s', mean(c1_fg), std(c1_fg), fname), ...
    'Interpreter','none', 'FontSize',10);
grid on; box on;
saveas(fig4, fullfile(out_dir, 'c1_distribution.png')); close(fig4);

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

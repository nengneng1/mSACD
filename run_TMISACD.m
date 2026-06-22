% =========================================================================
% run_TMISACD.m — TMI-SACD 双色超分辨完整管线
% =========================================================================
% 每个文件的结果独立存放在: result/{exp_name}/{name}/
%
% 对应原始管线: TMISACD_step1 → step2 → TMI → step3 → step4 → step5 → step6
% =========================================================================

% clc; clear; close all;

addpath('TMI_script/TMI_script/TMISACDF');
addpath('TMI_script/TMI_script/sparse/F');
addpath('TMI_script/TMI_script/pearson');
addpath('TMI_script/TMI_script/config');

%% ====== 配置（从本机 tmisacd_config.m 加载，不提交）======
% 复制 tmisacd_config_template.m → tmisacd_config.m 并填入本机路径。
tmisacd_config;
% 注入: base_dir, file_pattern, name_pattern, idx_list,
%       exp_name, confname, rsFPs, frame, skip

output_root = 'result';

%% ========================================================================
%  4. 算法参数
%  ========================================================================

% --- 预处理 ---
sigma       = [1, 1, 0];   % 高斯去噪 [sx, sy, sz], sz=0 为 2D 去噪
RL1_on      = 0;           % TMI 输入是否用 RL 解卷积
FWHM_TMI    = 3;           % TMI 阶段 PSF 半高宽 (像素)
iterTMI     = 10;          % TMI 阶段 RL 迭代次数

% --- TMI 拆分 ---
gamma_on    = 1;           % 是否 gamma 校正衰减曲线
gamma_1     = 1;           % gamma 值
bg_thresh   = 50/65535;    % 背景阈值 (低于此值跳过拟合)
mask_thresh = 1000/65535;  % 掩膜小值过滤阈值 (≈0.015, 设0可完全关闭)
FWHM4gauss  = [1, 1];      % mask 高斯平滑 [ch1, ch2]

% --- 超分辨 ---
finter          = 2;       % Fourier 插值倍数
paddingfactor   = 10;      % 边缘填充因子

% --- SACD RL 解卷积 ---
FWHM_sacd   = [3, 2.5];    % [第一轮, 第二轮] PSF 半高宽
iter_sacd   = [1, 7];      % [第一轮, 第二轮] 迭代次数

% --- SACD gamma 扫描 ---
gamma_scan  = [0.7, 1];

% --- 最终稀疏解卷积 (per channel) ---
%   mu: 保真度(越小越平滑), l1: 稀疏度, backg: 背景滤除(0=不滤)
sp_ch1_mu = 100;   sp_ch1_sigmat = 0;  sp_ch1_l1 = 0.01;  sp_ch1_iter = 100;  sp_ch1_backg = 0;
sp_ch2_mu = 100;   sp_ch2_sigmat = 0;  sp_ch2_l1 = 1;     sp_ch2_iter = 100;  sp_ch2_backg = 0;

FWHM_post   = 2.5;         % 最终 RL PSF 半高宽
iter_post   = 3;           % 最终 RL 迭代次数
order       = 2;           % PSF 阶数

%% ====== 初始化 ======

fprintf('共 %d 个文件待处理\n', numel(idx_list));
for i = 1:numel(idx_list)
    fprintf('  [%03d] %s\n', idx_list(i), sprintf(file_pattern, idx_list(i)));
end

load([confname, '.mat'], 'ch1', 'ch2');

%% ========================================================================
%  6. 逐文件处理
%  ========================================================================

for idx = idx_list

    fname     = sprintf(name_pattern, idx);
    imgname   = [base_dir, sprintf(file_pattern, idx)];

    % ---- 当前文件的输出目录 ----
    imgfolder = fullfile(output_root, exp_name, fname);
    if ~exist(imgfolder, 'dir'), mkdir(imgfolder); end

    fprintf('\n========== [%03d] %s ==========\n', idx, fname);
    fprintf('  输出: %s\n', imgfolder);

    %% --- 6.1 加载数据 ---
    data = double(imreadstack(imgname));
    n_total_frames = size(data, 3);
    n_SR_frames = floor(n_total_frames / skip);
    fprintf('  总帧数: %d, SR 组数: %d (skip=%d)\n', n_total_frames, n_SR_frames, skip);

    % 累加器: 每个 gamma 一组 {SACD1, SACD2, SOFI1, SOFI2}
    accum = cell(length(gamma_scan), 4);

    %% ====================================================================
    %  6.2 逐 SR 帧处理
    %  ====================================================================

    for f = 1:n_SR_frames
        fprintf('  SR 组 %d/%d\n', f, n_SR_frames);

        % ---- 提取帧 + 填充 + 归一化 ----
        imgstack = data(:, :, (f-1)*skip+1 : f*skip);
        imgstack = padarray(imgstack, [paddingfactor, paddingfactor, 0], 'symmetric');
        imgstack = imgstack(:, :, 1:frame);
        imgstack = imgstack ./ max(imgstack(:));
        imwritestack(uint16(imgstack .* 65535), fullfile(imgfolder, sprintf('%s_SR%d_raw.tif', fname, f)));

        % ---- 高斯去噪 ----
        if sigma(3) ~= 0
            imgstack_gauss = imgaussfilt3(imgstack, sigma);
        else
            imgstack_gauss = imgstack;
        end

        % ---- RL 解卷积 (TMI 预处理) ----
        IpsfTMI = generate_psfv0(FWHM_TMI);
        stack4TMI = zeros(size(imgstack_gauss));
        for i = 1:frame
            stack4TMI(:, :, i) = deconvlucy(imgstack_gauss(:, :, i), IpsfTMI, iterTMI);
        end

        if RL1_on == 1
            stack1 = stack4TMI;
        else
            stack1 = imgstack_gauss;
        end
        clear imgstack_gauss stack4TMI;

        % ---- Gamma 校正 ----
        if gamma_on == 1
            stack1_g = stack1 .^ gamma_1;
            ch1_g = (ch1 .^ gamma_1) ./ max(ch1 .^ gamma_1);
            ch2_g = (ch2 .^ gamma_1) ./ max(ch2 .^ gamma_1);
        else
            stack1_g = stack1; ch1_g = ch1; ch2_g = ch2;
        end

        % Reference curve separability check
        r1 = ch1_g(1:frame); r2 = ch2_g(1:frame);
        if isrow(r1), r1 = r1'; end
        if isrow(r2), r2 = r2'; end
        R_ref = [r1, r2];
        c_ref = corrcoef(R_ref);
        fprintf('  Reference: corr(ch1,ch2)=%.4f, cond([ch1,ch2])=%.2f\n', ...
            c_ref(1,2), cond(R_ref));

        % ---- PCA 去噪 (TMISACD_step2) ----
        stack1_norm = double(stack1_g ./ max(stack1_g(:)));
        filt_no_PC = 3;
        NoiseCorrection = sqrt(mean(stack1_norm, [1, 2]) + eps);
        temp_norm = stack1_norm ./ NoiseCorrection;
        temp_norm(temp_norm < 0) = 0;
        vector_conv = reshape(temp_norm, [size(temp_norm,1)*size(temp_norm,2), size(temp_norm,3)]);
        [flimcoeff, scores, ~] = pca(vector_conv);
        filt_vector_PCA = (scores(:, 1:filt_no_PC) * flimcoeff(:, 1:filt_no_PC)') + mean(vector_conv, 1);
        filt_img_PCA = reshape(filt_vector_PCA, [size(temp_norm,1), size(temp_norm,2), size(temp_norm,3)]);
        stack1_norm = (filt_img_PCA .* NoiseCorrection(:, :, :));
        stack1_norm = stack1_norm ./ max(stack1_norm(:));

        % 保存 PCA 去噪后的序列
        imwritestack(uint16(stack1_norm .* 65535), fullfile(imgfolder, sprintf('%s_SR%d_PCA_denoised.tif', fname, f)));

        % ================================================================
        %  TMI 解混: 同时运行旧(fmincon)和新(NNLS)方法, 仅用于对比
        %  后续 SACD 主流程使用 comp_raw_new (NNLS)
        % ================================================================
        [H, W, ~] = size(stack1_norm);
        cp = 1 + paddingfactor;

        % --- 旧方法: fmincon sum-to-one 比例解混 (仅对比, 不入后续流程) ---
        comp_raw_old = unmixTMI_FMINCON_old(stack1_norm, ch1_g, ch2_g, frame, bg_thresh, rsFPs);

        TMIch1_old = comp_raw_old(:, :, 1) .* imgstack(:, :, 1);
        TMIch2_old = comp_raw_old(:, :, 2) .* imgstack(:, :, 1);
        TMIch1_old = safe_norm01(TMIch1_old);
        TMIch2_old = safe_norm01(TMIch2_old);

        imwritestack(uint16(TMIch1_old(cp:end-paddingfactor, cp:end-paddingfactor) .* 65535), ...
            fullfile(imgfolder, sprintf('%s_SR%d_TMI_OLD_fmincon_ch1.tif', fname, f)));
        imwritestack(uint16(TMIch2_old(cp:end-paddingfactor, cp:end-paddingfactor) .* 65535), ...
            fullfile(imgfolder, sprintf('%s_SR%d_TMI_OLD_fmincon_ch2.tif', fname, f)));
        imwritestack(uint16(comp_raw_old(cp:end-paddingfactor, cp:end-paddingfactor, 1) .* 65535), ...
            fullfile(imgfolder, sprintf('%s_SR%d_ratio_OLD_fmincon_ch1.tif', fname, f)));
        imwritestack(uint16(comp_raw_old(cp:end-paddingfactor, cp:end-paddingfactor, 2) .* 65535), ...
            fullfile(imgfolder, sprintf('%s_SR%d_ratio_OLD_fmincon_ch2.tif', fname, f)));

        % --- 新方法: weighted NNLS 非负幅度解混 (后续 SACD 主流程) ---
        [amp_raw_new, comp_raw_new] = unmixTMI_NNLS(stack1_norm, ch1_g, ch2_g, frame, bg_thresh, rsFPs);
        clear stack1 stack1_g stack1_norm;

        TMIch1_new = safe_norm01(amp_raw_new(:, :, 1));
        TMIch2_new = safe_norm01(amp_raw_new(:, :, 2));

        imwritestack(uint16(TMIch1_new(cp:end-paddingfactor, cp:end-paddingfactor) .* 65535), ...
            fullfile(imgfolder, sprintf('%s_SR%d_TMI_NEW_NNLS_noB_amp_ch1.tif', fname, f)));
        imwritestack(uint16(TMIch2_new(cp:end-paddingfactor, cp:end-paddingfactor) .* 65535), ...
            fullfile(imgfolder, sprintf('%s_SR%d_TMI_NEW_NNLS_noB_amp_ch2.tif', fname, f)));
        imwritestack(uint16(amp_raw_new(cp:end-paddingfactor, cp:end-paddingfactor, 1) .* 65535), ...
            fullfile(imgfolder, sprintf('%s_SR%d_amp_NEW_NNLS_noB_ch1.tif', fname, f)));
        imwritestack(uint16(amp_raw_new(cp:end-paddingfactor, cp:end-paddingfactor, 2) .* 65535), ...
            fullfile(imgfolder, sprintf('%s_SR%d_amp_NEW_NNLS_noB_ch2.tif', fname, f)));
        imwritestack(uint16(comp_raw_new(cp:end-paddingfactor, cp:end-paddingfactor, 1) .* 65535), ...
            fullfile(imgfolder, sprintf('%s_SR%d_ratio_NEW_NNLS_noB_ch1.tif', fname, f)));
        imwritestack(uint16(comp_raw_new(cp:end-paddingfactor, cp:end-paddingfactor, 2) .* 65535), ...
            fullfile(imgfolder, sprintf('%s_SR%d_ratio_NEW_NNLS_noB_ch2.tif', fname, f)));

        % ---- 掩膜后处理 (TMISACD_step3, 使用 comp_raw_new) ----
        comp_img = min(max(comp_raw_new, 0), 1);
        comp_img(comp_img < mask_thresh) = 0;
        for ch = 1:rsFPs
            if FWHM4gauss(ch) ~= 0
                comp_img(:, :, ch) = imgaussfilt(comp_img(:, :, ch), FWHM4gauss(ch));
            end
        end
        comp_resized = imresize3(comp_img, ...
            [finter * H, finter * W, rsFPs], 'Method', 'cubic');
        imwritestack(uint16(comp_resized .* 65535), ...
            fullfile(imgfolder, sprintf('%s_SR%d_mask.tif', fname, f)));

        % ---- SACD RL 解卷积 (TMISACD_step4) ----
        stack4SACD = zeros(size(imgstack));
        Ipsf_sacd = generate_psfv0(FWHM_sacd(1));
        for i = 1:frame
            stack4SACD(:, :, i) = deconvlucy(imgstack(:, :, i), Ipsf_sacd, iter_sacd(1));
        end
        Ipsf_sacd2 = generate_psfv0(FWHM_sacd(2));
        for i = 1:frame
            stack4SACD(:, :, i) = deconvlucy(stack4SACD(:, :, i), Ipsf_sacd2, iter_sacd(2));
        end
        % 保存 SACD RL 解卷积后的序列 (原始 TMISACD_step4 输出)
        imwritestack(uint16(percennorm(stack4SACD, 0, 100) .* 65535), ...
            fullfile(imgfolder, sprintf('%s_SR%d_raw_RL.tif', fname, f)));

        sacd_resized = abs(fourierInterpolation(stack4SACD, [finter, finter, 1], 'lateral'));
        sacd_resized(sacd_resized < 0) = 0;   % 清除插值振铃负值
        sacd_resized = sacd_resized ./ max(sacd_resized(:));

        % ---- 3D 衰减曲线 (TMISACD_step5) ----
        ch1_ = zeros(1, 1, frame);
        ch2_ = zeros(1, 1, frame);
        for nn = 1:frame
            ch1_(1, 1, nn) = ch1(nn);
            ch2_(1, 1, nn) = ch2(nn);
        end

        % ---- SACD Cumulant + 最终重建 (TMISACD_step5 + step6) ----
        for ig = 1:length(gamma_scan)
            ga = gamma_scan(ig);

            maskimg1 = real(comp_resized(:, :, 1) .* ch1_ .* sacd_resized.^ga).^(1/ga);
            maskimg2 = real(comp_resized(:, :, 2) .* ch2_ .* sacd_resized.^ga).^(1/ga);
            maskimg1 = maskimg1 ./ max(maskimg1(:));
            maskimg2 = maskimg2 ./ max(maskimg2(:));
            for z = 1:size(maskimg1, 3)
                maskimg1z = maskimg1(:, :, z);
                maskimg2z = maskimg2(:, :, z);
                maskimg1(:, :, z) = maskimg1(:, :, z) ./ max(maskimg1z(:));
                maskimg2(:, :, z) = maskimg2(:, :, z) ./ max(maskimg2z(:));
            end

            stacksub1 = abs(maskimg1);   % mul=0
            stacksub2 = abs(maskimg2);
            stacksub1 = stacksub1(:, :, 1:20);
            stacksub2 = stacksub2(:, :, 1:20);

            % 保存细胞器序列中间结果
            tag_f = sprintf('SR%d_gamma%.1f', f, ga);
            imwritestack(uint16(maskimg1 .* 65535), fullfile(imgfolder, sprintf('%s_maskimg_ch1_%s.tif', fname, tag_f)));
            imwritestack(uint16(maskimg2 .* 65535), fullfile(imgfolder, sprintf('%s_maskimg_ch2_%s.tif', fname, tag_f)));
            imwritestack(uint16(stacksub1 .* 65535), fullfile(imgfolder, sprintf('%s_fluct_ch1_%s.tif', fname, tag_f)));
            imwritestack(uint16(stacksub2 .* 65535), fullfile(imgfolder, sprintf('%s_fluct_ch2_%s.tif', fname, tag_f)));

            % SOFI cumulant
            Nx = size(stacksub1, 1); Ny = size(stacksub1, 2);
            cum1 = zeros(Nx, Ny); cum2 = zeros(Nx, Ny);
            cum1(2:Nx-1, 2:Ny-1) = (mean(stacksub1(1:Nx-2, 2:Ny-1, :) .* stacksub1(3:Nx, 2:Ny-1, :), 3) ...
                + mean(stacksub1(2:Nx-1, 1:Ny-2, :) .* stacksub1(2:Nx-1, 3:Ny, :), 3)) ./ 2;
            cum2(2:Nx-1, 2:Ny-1) = (mean(stacksub2(1:Nx-2, 2:Ny-1, :) .* stacksub2(3:Nx, 2:Ny-1, :), 3) ...
                + mean(stacksub2(2:Nx-1, 1:Ny-2, :) .* stacksub2(2:Nx-1, 3:Ny, :), 3)) ./ 2;

            % 稀疏解卷积 + RL (TMISACD_step6)
            Ipsf2_1 = generate_psfv0(FWHM_post * finter);
            sparse1 = sparse_main(cum1.^0.5, sp_ch1_mu, sp_ch1_sigmat, sp_ch1_l1, sp_ch1_iter, sp_ch1_backg).^2;
            SACD1 = double(abs(deconvlucy(sparse1, Ipsf2_1.^order, iter_post)));
            SACD1 = double(SACD1 ./ max(SACD1(:)));
            SACD1 = SACD1.^0.5;
            SACD1 = SACD1 ./ max(SACD1(:));

            Ipsf2_2 = generate_psfv0(FWHM_post * finter);
            sparse2 = sparse_main(cum2.^0.5, sp_ch2_mu, sp_ch2_sigmat, sp_ch2_l1, sp_ch2_iter, sp_ch2_backg).^2;
            SACD2 = double(abs(deconvlucy(sparse2, Ipsf2_2.^order, iter_post)));
            SACD2 = double(SACD2 ./ max(SACD2(:)));
            SACD2 = SACD2.^0.5;
            SACD2 = SACD2 ./ max(SACD2(:));

            % 累加
            if isempty(accum{ig, 1})
                accum{ig, 1} = SACD1; accum{ig, 2} = SACD2;
                accum{ig, 3} = cum1;   accum{ig, 4} = cum2;
            else
                accum{ig, 1}(:, :, f) = SACD1; accum{ig, 2}(:, :, f) = SACD2;
                accum{ig, 3}(:, :, f) = cum1;   accum{ig, 4}(:, :, f) = cum2;
            end
        end

        clear imgstack stack4SACD sacd_resized comp_raw_new comp_raw_old amp_raw_new;
    end  % f loop

    %% --- 6.3 保存当前文件的最终输出 ---
    c1 = paddingfactor * finter + 1;
    crop = @(img) img(c1:end-paddingfactor*finter, c1:end-paddingfactor*finter, :);

    for ig = 1:length(gamma_scan)
        ga = gamma_scan(ig);
        tag = sprintf('gamma%.1f', ga);
        imwritestack(uint16(crop(accum{ig, 1}) .* 65535), fullfile(imgfolder, sprintf('%s_TMISACD_ch1_%s.tif', fname, tag)));
        imwritestack(uint16(crop(accum{ig, 2}) .* 65535), fullfile(imgfolder, sprintf('%s_TMISACD_ch2_%s.tif', fname, tag)));
        imwritestack(uint16(crop(accum{ig, 3}) .* 65535), fullfile(imgfolder, sprintf('%s_SOFI_ch1_%s.tif', fname, tag)));
        imwritestack(uint16(crop(accum{ig, 4}) .* 65535), fullfile(imgfolder, sprintf('%s_SOFI_ch2_%s.tif', fname, tag)));
    end

    fprintf('========== 完成 %s ==========\n', fname);
end

fprintf('\n===== 全部 %d 个文件处理完成 =====\n', numel(idx_list));

% =========================================================================
function [amp_raw, comp_raw] = unmixTMI_NNLS(stack1_norm, ch1_g, ch2_g, frame, bg_thresh, rsFPs)
% Weighted NNLS unmixing (no baseline):  curve(t) ≈ a1*ch1(t) + a2*ch2(t)
%   a1 >= 0, a2 >= 0
%
% Weighting: w = 1/sqrt(curve+1e-4)  (Poisson-like, 避免亮帧主导 loss)
%
% Input:
%   stack1_norm - H×W×F  normalized intensity stack
%   ch1_g, ch2_g - reference decay curves (column or row vectors)
%   frame       - number of time frames
%   bg_thresh   - background threshold
%   rsFPs       - number of fluorophores (2)
%
% Output:
%   amp_raw  - H×W×rsFPs  [a1, a2] fitted absolute amplitudes
%   comp_raw - H×W×rsFPs  [c1, c2] = a_i / sum(a_i+eps), for mask downstream

[H, W, F] = size(stack1_norm);

% Design matrix: [ch1, ch2]  (frame × 2, no baseline)
c1 = ch1_g(1:frame); if isrow(c1), c1 = c1'; end
c2 = ch2_g(1:frame); if isrow(c2), c2 = c2'; end
R = [c1, c2];

amp_raw = zeros(H, W, rsFPs);
comp_raw = zeros(H, W, rsFPs);

parfor ii = 1:H
    amp_row = zeros(W, rsFPs);
    comp_row = zeros(W, rsFPs);
    for jj = 1:W
        curve = squeeze(double(stack1_norm(ii, jj, :)));
        if isrow(curve), curve = curve'; end

        curve(curve < 0) = 0;  % PCA 去噪可能引入微小负值, 截断

        if max(curve) < bg_thresh
            amp_row(jj, :) = 0;
            comp_row(jj, :) = 0;
        else
            % Poisson-like weighting
            w = 1 ./ sqrt(curve + 1e-4);
            Rw = R .* w;
            yw = curve .* w;

            % Non-negative least squares
            coef = lsqnonneg(Rw, yw);

            % a1, a2 (coef(3) = baseline, discarded)
            amp_row(jj, :) = coef(1:rsFPs)';

            % Proportions for mask
            s = sum(coef(1:rsFPs)) + eps;
            comp_row(jj, :) = coef(1:rsFPs)' / s;
        end
    end
    amp_raw(ii, :, :) = amp_row;
    comp_raw(ii, :, :) = comp_row;
end
end

% =========================================================================
function comp_raw_old = unmixTMI_FMINCON_old(stack1_norm, ch1_g, ch2_g, frame, bg_thresh, rsFPs)
% Old TMI unmixing via fmincon:  curve ≈ c1*ch1 + c2*ch2
%   Constraints: c1 + c2 = 1,  0 <= c1, c2 <= 1
%
% 此函数仅用于与 NNLS 方法对比, 不参与后续 SACD 流程.

[H, W, F] = size(stack1_norm); %#ok<NASGU>

c1_ref = ch1_g(1:frame); if isrow(c1_ref), c1_ref = c1_ref'; end
c2_ref = ch2_g(1:frame); if isrow(c2_ref), c2_ref = c2_ref'; end

comp_raw_old = zeros(H, W, rsFPs);
options = optimoptions('fmincon', 'Display', 'none');

parfor ii = 1:H
    comp_row = zeros(W, rsFPs);
    for jj = 1:W
        curve = squeeze(double(stack1_norm(ii, jj, :)));
        if isrow(curve), curve = curve'; end
        curve(curve < 0) = 0;

        if max(curve) < bg_thresh
            comp_row(jj, :) = 0;
        else
            [pv, ~, ~, ~] = fmincon( ...
                @(c) signalseparation2fmincon(c1_ref, c2_ref, curve, c), ...
                [0.5, 0.5], [], [], ones(1, rsFPs), 1, ...
                zeros(1, rsFPs), ones(1, rsFPs), [], options);
            comp_row(jj, :) = pv(1:rsFPs);
        end
    end
    comp_raw_old(ii, :, :) = comp_row;
end
end

% =========================================================================
function out = safe_norm01(img)
% Safe normalization to [0,1], returns 0 if all values are 0.
m = max(img(:));
if m > 0
    out = img ./ m;
else
    out = img;
end
end

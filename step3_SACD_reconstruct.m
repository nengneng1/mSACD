% =========================================================================
% step3_SACD_reconstruct.m — SACD 超分辨重建
% =========================================================================
% 输入 (workspace): imgstack, comp_resized, frame, paddingfactor, finter,
%                   FILE_OUT_DIR, fname, skip, data_all, n_SR_frames
% 输出 (tif):       {FILE_OUT_DIR}/SACD/{fname}_TMISACD_chX_gammaX.X.tif 等
% =========================================================================

%% ====== 路径 ======
STEP_NAME = 'SACD';
addpath(genpath('F'));

%% ====== 参数 ======
% SACD RL (TMISACD_step4)
FWHM_sacd  = [3];
iter_sacd  = [7];

% SACD cumulant (TMISACD_step5)
order      = 2;
gamma_scan = [0.5,0.7,1];

% 稀疏解卷积 (TMISACD_step6)
sp_ch1_mu = 500;  sp_ch1_sigmat = 0;  sp_ch1_l1 = 1;  sp_ch1_iter = 100;  sp_ch1_backg = 0;
sp_ch2_mu = 500;  sp_ch2_sigmat = 0;  sp_ch2_l1 = 1;  sp_ch2_iter = 100;  sp_ch2_backg = 0;
FWHM_post  = 3;
iter_post  = 3;

%% ====== 输出目录 ======
out_dir = fullfile(FILE_OUT_DIR, STEP_NAME);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ====== 衰减曲线 ======
load([confname, '.mat'], 'ch1', 'ch2');

%% ====== 逐 SR 帧处理 ======
accum = cell(length(gamma_scan), 4);
paddingfactor_step3 = paddingfactor;  % 避免与其他变量冲突

for f = 1:n_SR_frames
    fprintf('  SR 组 %d/%d\n', f, n_SR_frames);

    % --- 提取帧 ---
    imgstack_f = data_all(:, :, (f-1)*skip+1 : f*skip);
    imgstack_f = padarray(imgstack_f, [paddingfactor_step3, paddingfactor_step3, 0], 'symmetric');
    imgstack_f = imgstack_f(:, :, 1:frame);
    imgstack_f = imgstack_f ./ max(imgstack_f(:));

    % --- SACD RL 解卷积 (TMISACD_step4) ---
    stack4SACD = zeros(size(imgstack_f));
    Ipsf1 = generate_psfv0(FWHM_sacd(1));
    for i = 1:frame
        stack4SACD(:, :, i) = deconvlucy(imgstack_f(:, :, i), Ipsf1, iter_sacd(1));
    end


    % 保存 SACD RL 解卷积后的序列 (原始 TMISACD_step4 输出)
    imwritestack(uint16(percennorm(stack4SACD, 0, 100) .* 65535), ...
        fullfile(out_dir, sprintf('SR%d_raw_RL.tif', f)));

    sacd_resized = abs(fourierInterpolation(stack4SACD, [finter, finter, 1], 'lateral'));
    sacd_resized(sacd_resized < 0) = 0;
    sacd_resized = sacd_resized ./ max(sacd_resized(:));
    clear stack4SACD;

    % --- 3D 衰减曲线 ---
    ch1_ = zeros(1, 1, frame); ch2_ = zeros(1, 1, frame);
    for nn = 1:frame
        ch1_(1, 1, nn) = ch1(nn); ch2_(1, 1, nn) = ch2(nn);
    end

    % --- 多 gamma SACD (TMISACD_step5 + step6) ---
    for ig = 1:length(gamma_scan)
        ga = gamma_scan(ig);

        maskimg1 = real(comp_resized(:, :, 1) .* ch1_ .* sacd_resized.^ga).^(1/ga);
        maskimg2 = real(comp_resized(:, :, 2) .* ch2_ .* sacd_resized.^ga).^(1/ga);
        maskimg1 = maskimg1 ./ max(maskimg1(:));
        maskimg2 = maskimg2 ./ max(maskimg2(:));

        stacksub1 = abs(maskimg1); stacksub1 = stacksub1(:, :, 1:20);
        stacksub2 = abs(maskimg2); stacksub2 = stacksub2(:, :, 1:20);

        % 保存中间序列 (mask × 数据 = 单个细胞器的时序信号)
        tag_f = sprintf('SR%d_gamma%.1f', f, ga);
        imwritestack(uint16(maskimg1 .* 65535), fullfile(out_dir, sprintf('maskimg_ch1_%s.tif', tag_f)));
        imwritestack(uint16(maskimg2 .* 65535), fullfile(out_dir, sprintf('maskimg_ch2_%s.tif', tag_f)));
        imwritestack(uint16(stacksub1 .* 65535), fullfile(out_dir, sprintf('fluct_ch1_%s.tif', tag_f)));
        imwritestack(uint16(stacksub2 .* 65535), fullfile(out_dir, sprintf('fluct_ch2_%s.tif', tag_f)));

        Nx = size(stacksub1, 1); Ny = size(stacksub1, 2);
        cum1 = zeros(Nx, Ny); cum2 = zeros(Nx, Ny);
        cum1(2:Nx-1, 2:Ny-1) = (mean(stacksub1(1:Nx-2, 2:Ny-1, :) .* stacksub1(3:Nx, 2:Ny-1, :), 3) ...
            + mean(stacksub1(2:Nx-1, 1:Ny-2, :) .* stacksub1(2:Nx-1, 3:Ny, :), 3)) ./ 2;
        cum2(2:Nx-1, 2:Ny-1) = (mean(stacksub2(1:Nx-2, 2:Ny-1, :) .* stacksub2(3:Nx, 2:Ny-1, :), 3) ...
            + mean(stacksub2(2:Nx-1, 1:Ny-2, :) .* stacksub2(2:Nx-1, 3:Ny, :), 3)) ./ 2;

        % ch1 稀疏 + RL
        sp1 = sparse_main(cum1.^0.5, sp_ch1_mu, sp_ch1_sigmat, sp_ch1_l1, sp_ch1_iter, sp_ch1_backg).^2;
        Ipsf_p = generate_psfv0(FWHM_post * finter);
        SACD1 = double(abs(deconvlucy(sp1, Ipsf_p.^order, iter_post)));
        SACD1 = double(SACD1 ./ max(SACD1(:)));
        SACD1 = SACD1.^0.5;
        SACD1 = SACD1 ./ max(SACD1(:));

        % ch2 稀疏 + RL
        sp2 = sparse_main(cum2.^0.5, sp_ch2_mu, sp_ch2_sigmat, sp_ch2_l1, sp_ch2_iter, sp_ch2_backg).^2;
        SACD2 = double(abs(deconvlucy(sp2, Ipsf_p.^order, iter_post)));
        SACD2 = double(SACD2 ./ max(SACD2(:)));
        SACD2 = SACD2.^0.5;
        SACD2 = SACD2 ./ max(SACD2(:));

        if isempty(accum{ig, 1})
            accum{ig, 1} = SACD1; accum{ig, 2} = SACD2;
            accum{ig, 3} = cum1;   accum{ig, 4} = cum2;
        else
            accum{ig, 1}(:, :, f) = SACD1; accum{ig, 2}(:, :, f) = SACD2;
            accum{ig, 3}(:, :, f) = cum1;   accum{ig, 4}(:, :, f) = cum2;
        end
    end
    clear imgstack_f sacd_resized;
end

%% ====== 保存 ======
c1 = paddingfactor_step3 * finter + 1;
crop = @(x) x(c1:end-paddingfactor_step3*finter, c1:end-paddingfactor_step3*finter, :);

for ig = 1:length(gamma_scan)
    ga = gamma_scan(ig);
    tag = sprintf('gamma%.1f', ga);
    imwritestack(uint16(crop(accum{ig, 1}) .* 65535), fullfile(out_dir, sprintf('TMISACD_ch1_%s.tif', tag)));
    imwritestack(uint16(crop(accum{ig, 2}) .* 65535), fullfile(out_dir, sprintf('TMISACD_ch2_%s.tif', tag)));
    imwritestack(uint16(crop(accum{ig, 3}) .* 65535), fullfile(out_dir, sprintf('SOFI_ch1_%s.tif', tag)));
    imwritestack(uint16(crop(accum{ig, 4}) .* 65535), fullfile(out_dir, sprintf('SOFI_ch2_%s.tif', tag)));
end

fprintf('  [step3] SACD重建完成 → %s\n', out_dir);

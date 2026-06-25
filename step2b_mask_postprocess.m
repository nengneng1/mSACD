% =========================================================================
% step2b_mask_postprocess.m — 掩膜后处理 (阈值 + 高斯平滑 + 上采样)
% =========================================================================
% 输入 (workspace): comp_raw, H, W, FILE_OUT_DIR, fname
% 输出 (workspace): comp_img, comp_resized, finter
% 输出 (tif):       {FILE_OUT_DIR}/mask/{fname}_mask.tif
% =========================================================================

%% ====== 路径 ======
STEP_NAME = 'mask';
addpath(genpath('F'));

%% ====== 参数 ======
rsFPs       = 2;
mask_thresh = 0/65535;
FWHM4gauss  = [1, 1];
finter      = 2;

%% ====== 输出目录 (按解混方法标签区分, 避免互相覆盖) ======
if ~exist('UNMIX_TAG', 'var') || isempty(UNMIX_TAG)
    out_dir = fullfile(FILE_OUT_DIR, STEP_NAME);
else
    out_dir = fullfile(FILE_OUT_DIR, [STEP_NAME '_' UNMIX_TAG]);
end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ====== 掩膜后处理 (TMISACD_step3) ======
comp_img = min(max(comp_raw, 0), 1);
comp_img(comp_img < mask_thresh) = 0;

for ch = 1:rsFPs
    if FWHM4gauss(ch) ~= 0
        comp_img(:, :, ch) = imgaussfilt(comp_img(:, :, ch), FWHM4gauss(ch));
    end
end

comp_resized = imresize3(comp_img, ...
    [finter * H, finter * W, rsFPs], 'Method', 'cubic');

imwritestack(uint16(comp_resized .* 65535), fullfile(out_dir, 'mask.tif'));

fprintf('  [step2b] mask后处理完成 → %s\n', out_dir);

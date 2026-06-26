% =========================================================================
% step2_TMI_unmix.m — PCA去噪 + TMI 通道拆分 (fmincon sum-to-one)
% =========================================================================
% 输入 (workspace): stack1_for_TMI, imgstack, frame, paddingfactor,
%                   FILE_OUT_DIR, fname
% 输出 (workspace): amp_raw, comp_raw, TMIch1, TMIch2, H, W
% 输出 (tif):       {FILE_OUT_DIR}/TMI_unmix/
%
% 解混模型: curve ≈ c1*ch1 + c2*ch2,  c1+c2=1,  0<=c1,c2<=1
% =========================================================================

%% ====== 路径 ======
STEP_NAME = 'TMI_unmix';
UNMIX_TAG = 'unmix';     % 方法标签 → 供 step2b/step3 区分输出目录
addpath(genpath('F'));

%% ====== 参数 ======
rsFPs     = 2;
bg_thresh = 50/65535;
pca_on    = 0;          % PCA 去噪开关 (0=关闭, 用于对比)

%% ====== 输出目录 ======
out_dir = fullfile(FILE_OUT_DIR, STEP_NAME);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ====== 衰减曲线 ======
load([confname, '.mat'], 'ch1', 'ch2');

%% ====== 参考曲线 (无 gamma, 保持线性混合假设) ======
stack1_g = stack1_for_TMI; ch1_g = ch1; ch2_g = ch2;

% Reference curve separability check
R_ref = [ch1_g(1:frame), ch2_g(1:frame)];
c_ref = corrcoef(R_ref);
fprintf('  Reference: corr(ch1,ch2)=%.4f, cond([ch1,ch2])=%.2f\n', ...
    c_ref(1,2), cond(R_ref));

%% ====== PCA 去噪 (TMISACD_step2; pca_on=0 时跳过) ======
stack1_norm = double(stack1_g ./ max(stack1_g(:)));
if pca_on == 1
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
    imwritestack(uint16(stack1_norm .* 65535), fullfile(out_dir, 'PCA_denoised.tif'));
else
    imwritestack(uint16(stack1_norm .* 65535), fullfile(out_dir, 'noPCA_input.tif'));
end

%% ====== TMI fmincon 解混 ======
[H, W, ~] = size(stack1_norm);
comp_raw = unmixTMI_FMINCON(stack1_norm, ch1_g, ch2_g, frame, bg_thresh, rsFPs);
amp_raw  = comp_raw;  % fmincon 输出为比例，此处 amp_raw = comp_raw 供后续步骤兼容

%% ====== TMI 输出通道 ======
TMIch1 = comp_raw(:, :, 1) .* stack1_for_TMI(:, :, 1);
TMIch2 = comp_raw(:, :, 2) .* stack1_for_TMI(:, :, 1);
TMIch1 = TMIch1 ./ max(TMIch1(:) + eps);
TMIch2 = TMIch2 ./ max(TMIch2(:) + eps);

cp = 1 + paddingfactor;
crop = @(x) x(cp:end-paddingfactor, cp:end-paddingfactor, :);

imwritestack(uint16(crop(TMIch1) .* 65535), fullfile(out_dir, 'TMIch1.tif'));
imwritestack(uint16(crop(TMIch2) .* 65535), fullfile(out_dir, 'TMIch2.tif'));

fprintf('  [step2] TMI (fmincon) 拆分完成 → %s\n', out_dir);
fprintf('    comp_raw ch1: [%.4f, %.4f]  ch2: [%.4f, %.4f]\n', ...
    min(comp_raw(:,:,1),[],'all'), max(comp_raw(:,:,1),[],'all'), ...
    min(comp_raw(:,:,2),[],'all'), max(comp_raw(:,:,2),[],'all'));

% =========================================================================
function comp_raw = unmixTMI_FMINCON(stack1_norm, ch1_g, ch2_g, frame, bg_thresh, rsFPs)
% TMI unmixing via fmincon:  curve ≈ c1*ch1 + c2*ch2
%   Constraints: c1 + c2 = 1,  0 <= c1, c2 <= 1

[H, W, ~] = size(stack1_norm);
stack1_norm = real(stack1_norm);

c1_ref = real(ch1_g(1:frame)); if isrow(c1_ref), c1_ref = c1_ref'; end
c2_ref = real(ch2_g(1:frame)); if isrow(c2_ref), c2_ref = c2_ref'; end

comp_raw = zeros(H, W, rsFPs);
options  = optimoptions('fmincon', 'Display', 'none');

parfor ii = 1:H
    comp_row = zeros(W, rsFPs);
    for jj = 1:W
        if stack1_norm(ii, jj, 1) < bg_thresh
            comp_row(jj, :) = 0;
        else
            curve = real(squeeze(double(stack1_norm(ii, jj, :))));
            if isrow(curve), curve = curve'; end
            pv = fmincon( ...
                @(c) signalseparation2fmincon(c1_ref, c2_ref, curve, c), ...
                [0.5, 0.5], [], [], ones(1, rsFPs), 1, ...
                zeros(1, rsFPs), ones(1, rsFPs), [], options);
            comp_row(jj, :) = pv;
        end
    end
    comp_raw(ii, :, :) = comp_row;
end
end

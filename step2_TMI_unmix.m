% =========================================================================
% step2_TMI_unmix.m — PCA去噪 + TMI 通道拆分 (weighted NNLS)
% =========================================================================
% 输入 (workspace): stack1_for_TMI, imgstack, frame, paddingfactor,
%                   FILE_OUT_DIR, fname
% 输出 (workspace): amp_raw, comp_raw, TMIch1, TMIch2, H, W
% 输出 (tif):       {FILE_OUT_DIR}/TMI_unmix/
%
% 解混模型: curve ≈ a1*ch1 + a2*ch2 + b,  a1,a2,b >= 0
% 不再强制 c1+c2=1, 弱通道不再被亮通道吞掉
% =========================================================================

%% ====== 路径 ======
STEP_NAME = 'TMI_unmix';
addpath(genpath('F'));

%% ====== 参数 ======
rsFPs     = 2;
gamma_on  = 1;
gamma_1   = 1;
bg_thresh = 50/65535;

%% ====== 输出目录 ======
out_dir = fullfile(FILE_OUT_DIR, STEP_NAME);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ====== 衰减曲线 ======
load([confname, '.mat'], 'ch1', 'ch2');

%% ====== Gamma 校正 ======
if gamma_on == 1
    stack1_g = stack1_for_TMI .^ gamma_1;
    ch1_g = (ch1 .^ gamma_1) ./ max(ch1 .^ gamma_1);
    ch2_g = (ch2 .^ gamma_1) ./ max(ch2 .^ gamma_1);
else
    stack1_g = stack1_for_TMI; ch1_g = ch1; ch2_g = ch2;
end

% Reference curve separability check
R_ref = [ch1_g(1:frame), ch2_g(1:frame)];
c_ref = corrcoef(R_ref);
fprintf('  Reference: corr(ch1,ch2)=%.4f, cond([ch1,ch2])=%.2f\n', ...
    c_ref(1,2), cond(R_ref));

%% ====== PCA 去噪 (TMISACD_step2) ======
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
imwritestack(uint16(stack1_norm .* 65535), fullfile(out_dir, sprintf('%s_PCA_denoised.tif', fname)));

%% ====== TMI weighted NNLS 解混 ======
[H, W, ~] = size(stack1_norm);
[amp_raw, comp_raw] = unmixTMI_NNLS(stack1_norm, ch1_g, ch2_g, frame, bg_thresh, rsFPs);

%% ====== TMI 输出通道 (使用 amp_raw, 各自归一化) ======
TMIch1 = amp_raw(:, :, 1);
TMIch2 = amp_raw(:, :, 2);
TMIch1 = TMIch1 ./ max(TMIch1(:));
TMIch2 = TMIch2 ./ max(TMIch2(:));

cp = 1 + paddingfactor;
crop = @(x) x(cp:end-paddingfactor, cp:end-paddingfactor, :);

imwritestack(uint16(crop(TMIch1) .* 65535), fullfile(out_dir, sprintf('%s_TMIch1.tif', fname)));
imwritestack(uint16(crop(TMIch2) .* 65535), fullfile(out_dir, sprintf('%s_TMIch2.tif', fname)));

fprintf('  [step2] TMI (NNLS) 拆分完成 → %s\n', out_dir);
fprintf('    amp_raw ch1: [%.4f, %.4f]  ch2: [%.4f, %.4f]\n', ...
    min(amp_raw(:,:,1),[],'all'), max(amp_raw(:,:,1),[],'all'), ...
    min(amp_raw(:,:,2),[],'all'), max(amp_raw(:,:,2),[],'all'));

% =========================================================================
function [amp_raw, comp_raw] = unmixTMI_NNLS(stack1_norm, ch1_g, ch2_g, frame, bg_thresh, rsFPs)
% Weighted NNLS unmixing:  curve(t) ≈ a1*ch1(t) + a2*ch2(t) + b
%   a1 >= 0, a2 >= 0, b >= 0
%
% Weighting: w = 1/sqrt(curve+1e-4)  (Poisson-like, 避免亮帧主导 loss)
%
% Input:
%   stack1_norm - H×W×F  normalized intensity stack
%   ch1_g, ch2_g - column vectors, length=frame, reference decay curves
%   frame       - number of time frames used
%   bg_thresh   - background threshold (skip pixels below this)
%   rsFPs       - number of fluorophores (2)
%
% Output:
%   amp_raw  - H×W×rsFPs  absolute fitted amplitudes [a1, a2]
%   comp_raw - H×W×rsFPs  normalized proportions a_i / (a1+a2+eps)

[H, W, F] = size(stack1_norm);
stack1_norm = real(stack1_norm);

% Design matrix: [ch1, ch2, baseline=1]  (frame × 3)
c1 = ch1_g(1:frame); if isrow(c1), c1 = c1'; end
c2 = ch2_g(1:frame); if isrow(c2), c2 = c2'; end
R = [c1, c2, ones(frame, 1)];

amp_raw = zeros(H, W, rsFPs);
comp_raw = zeros(H, W, rsFPs);

parfor ii = 1:H
    amp_row = zeros(W, rsFPs);
    comp_row = zeros(W, rsFPs);
    for jj = 1:W
        if stack1_norm(ii, jj, 1) < bg_thresh
            amp_row(jj, :) = 0;
            comp_row(jj, :) = 0;
        else
            curve = squeeze(double(stack1_norm(ii, jj, :)));
            if isrow(curve), curve = curve'; end

            % Poisson-like weighting
            w = 1 ./ sqrt(curve + 1e-4);
            Rw = R .* w;
            yw = curve .* w;

            % Non-negative least squares
            coef = lsqnonneg(Rw, yw);

            % a1, a2 (skip baseline b = coef(3))
            amp_row(jj, :) = coef(1:rsFPs)';

            % Proportions (for mask downstream)
            s = sum(coef(1:rsFPs)) + eps;
            comp_row(jj, :) = coef(1:rsFPs)' / s;
        end
    end
    amp_raw(ii, :, :) = amp_row;
    comp_raw(ii, :, :) = comp_row;
end
end

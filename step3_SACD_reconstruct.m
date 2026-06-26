% =========================================================================
% step3_SACD_reconstruct.m — SACD 超分辨重建
% =========================================================================
% 输入 (workspace): comp_img, stack1_for_TMI, frame, paddingfactor,
%                   FILE_OUT_DIR, fname, n_SR_frames, confname
%                   (stack1_for_TMI: PCA 之前的预处理时序栈, 含 RL 与否取决于 step1 开关)
%                   (comp_img: H×W×2, 与 stack1_for_TMI 同分辨率)
% 输出 (tif):       {FILE_OUT_DIR}/SACD_{UNMIX_TAG}/
% =========================================================================
% Pipeline (每个 SR 帧组):
%   1. 时变分离        → mask×衰减曲线权重 → sig1, sig2
%   2. Dark sectioning → 暗通道背景估计并减除
%   3. 第一次 RL       → 每通道独立 RL 解卷积
%   4. SOFI 2 阶累积量 → 邻像素互相关
%   5. sparse_main     → 稀疏约束 (cum^0.5 → 平方)
%   6. 3× Fourier 上采样
%   7. 第二次 RL       → PSF 缩放到上采样后像素尺度
% =========================================================================

%% ====== 路径 ======
STEP_NAME = 'SACD';
addpath(genpath('F'));

%% ====== 参数 ======

% --- 第一次 RL [ch1, ch2] ---
FWHM_sacd = [2.8,  2.8];   % PSF 半高宽 (像素)
iter_sacd = [7, 5];   % 迭代次数

% --- Dark sectioning 去背景 ---
% 光学参数 (两通道共用同一光学系统)
ds_NA         = 1.4;
ds_emwl       = 600;    % 发射波长 (nm)
ds_pixelsize  = 110;    % 像素尺寸 (nm/px)
ds_factor     = 1;      % PSF 超采样因子
ds_pad        = 1;      % 1=对称填充边缘, 0=零填充
ds_divide     = 0.7;    % Hi/Lo 分频点
% 去背景参数 [ch1, ch2] 分别调节
ds_background = [0,    0   ];  % 0=轻度, 1=强力
ds_thres      = [70,   70  ];  % 背景阈值 (0–255 归一化尺度)
ds_dark_blend = [1,    1   ];  % 混合比 (1=完全去背景, 0=不变)
ds_pad_size   = [15,   15  ];  % 边缘填充分母
ds_dark3d     = [false, false]; % true=时序均值估背景, false=逐帧独立

% --- SOFI ---
n_sofi_frames = 20;     % 用于累积量计算的帧数

% --- sparse_main ---
sp_ch1_mu = 2000;  sp_ch1_sigmat = 0;  sp_ch1_l1 = 0.01;  sp_ch1_iter = 100;  sp_ch1_backg = 0;
sp_ch2_mu = 2000;  sp_ch2_sigmat = 0;  sp_ch2_l1 = 0.01;  sp_ch2_iter = 100;  sp_ch2_backg = 0;

% --- 上采样 + 第二次 RL [ch1, ch2] ---
finter_sr  = 2;              % Fourier 上采样倍数 (两通道共用)
FWHM_post  = [2.7,2.7];   % 第二次 RL PSF FWHM (原始像素单位)
iter_post  = [5,    5 ];   % 迭代次数
order      = [2,    2  ];   % cumulant 阶数 (deconvlucy 用 PSF^order)

%% ====== 输出目录 ======
if ~exist('UNMIX_TAG', 'var') || isempty(UNMIX_TAG)
    out_dir = fullfile(FILE_OUT_DIR, STEP_NAME);
else
    out_dir = fullfile(FILE_OUT_DIR, [STEP_NAME '_' UNMIX_TAG]);
end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ====== 衰减曲线 ======
load([confname, '.mat'], 'ch1', 'ch2');
ch1_t = reshape(real(ch1(1:frame)), 1, 1, frame);   % 1×1×frame (广播用)
ch2_t = reshape(real(ch2(1:frame)), 1, 1, frame);

%% ====== Dark sectioning 基础参数结构 ======
ds_params.NA           = ds_NA;
ds_params.emwavelength = ds_emwl;
ds_params.pixelsize    = ds_pixelsize;
ds_params.factor       = ds_factor;
% Nx/Ny 在 darksection_stack 内按实际图像尺寸填写

%% ====== 逐 SR 帧组处理 ======
accum_SACD1 = [];  accum_SACD2 = [];
accum_SOFI1 = [];  accum_SOFI2 = [];

% for f = 1:n_SR_frames
    fprintf('  SR 帧 %d/%d\n', f, n_SR_frames);

    %% 1. 时序数据（用 PCA 之前的预处理栈 stack1_for_TMI）
    % stack1_for_TMI 来自 step1: 已 padding/归一化的前 frame 帧, 是否含 RL 取决于 RL1_on 开关。
    % 不用 PCA 后的数据(会破坏时变分离), 也不用 data_all(已在 step1 完成 padding/RL 选择)。
    imgstack_f = double(stack1_for_TMI(:, :, 1:frame));
    imgstack_f = imgstack_f ./ max(imgstack_f(:) + eps);

    [Hf, Wf, ~] = size(imgstack_f);

    %% 2. 时变分离 → sig1, sig2
    % comp_img (H×W×2) 来自 step2b, 与 padded imgstack_f 同分辨率
    % 若尺寸不一致则双三次插值匹配
    if size(comp_img, 1) ~= Hf || size(comp_img, 2) ~= Wf
        comp_f = imresize3(comp_img, [Hf, Wf, 2], 'Method', 'cubic');
    else
        comp_f = comp_img;
    end
    comp_f = max(comp_f, 0);

    % 时变权重分母: H×W×frame (广播)
    denom = comp_f(:,:,1) .* ch1_t + comp_f(:,:,2) .* ch2_t + eps;
    sig1  = imgstack_f .* (comp_f(:,:,1) .* ch1_t) ./ denom;
    sig2  = imgstack_f .* (comp_f(:,:,2) .* ch2_t) ./ denom;
    sig1  = sig1 ./ max(sig1(:) + eps);
    sig2  = sig2 ./ max(sig2(:) + eps);

    imwritestack(uint16(sig1 .* 65535), fullfile(out_dir, sprintf('SR%d_sig1_sep.tif', f)));
    imwritestack(uint16(sig2 .* 65535), fullfile(out_dir, sprintf('SR%d_sig2_sep.tif', f)));

    %% 3. Dark sectioning 去背景（解混前不再做, 统一在此处对分离后的两通道执行）
    sig1 = darksection_stack(sig1, ds_params, ds_background(1), ds_pad, ds_thres(1), ...
                              ds_divide, ds_dark_blend(1), ds_pad_size(1), ds_dark3d(1));
    sig2 = darksection_stack(sig2, ds_params, ds_background(2), ds_pad, ds_thres(2), ...
                              ds_divide, ds_dark_blend(2), ds_pad_size(2), ds_dark3d(2));
    sig1 = sig1 ./ max(sig1(:) + eps);
    sig2 = sig2 ./ max(sig2(:) + eps);

    imwritestack(uint16(sig1 .* 65535), fullfile(out_dir, sprintf('SR%d_sig1_bgsub.tif', f)));
    imwritestack(uint16(sig2 .* 65535), fullfile(out_dir, sprintf('SR%d_sig2_bgsub.tif', f)));

    %% 4. 第一次 RL 解卷积
    Ipsf_ch1  = generate_psfv0(FWHM_sacd(1));
    Ipsf_ch2  = generate_psfv0(FWHM_sacd(2));
    RL1       = zeros(Hf, Wf, frame);
    RL2       = zeros(Hf, Wf, frame);
    for i = 1:frame
        RL1(:,:,i) = deconvlucy(sig1(:,:,i), Ipsf_ch1, iter_sacd(1));
        RL2(:,:,i) = deconvlucy(sig2(:,:,i), Ipsf_ch2, iter_sacd(2));
    end
    RL1 = RL1 ./ max(RL1(:) + eps);
    RL2 = RL2 ./ max(RL2(:) + eps);

    imwritestack(uint16(RL1 .* 65535), fullfile(out_dir, sprintf('SR%d_RL1_ch1.tif', f)));
    imwritestack(uint16(RL2 .* 65535), fullfile(out_dir, sprintf('SR%d_RL1_ch2.tif', f)));

    %% 5. SOFI 2 阶累积量 (邻像素互相关)
    nf   = min(n_sofi_frames, frame);
    sub1 = RL1(:,:,1:nf);
    sub2 = RL2(:,:,1:nf);
sub1 = abs(sub1-0.8*mean(sub1,3));
sub2 = abs(sub2-0.8*mean(sub2,3));
    cum1 = zeros(Hf, Wf);
    cum2 = zeros(Hf, Wf);
    cum1(2:Hf-1, 2:Wf-1) = (mean(sub1(1:Hf-2, 2:Wf-1,:) .* sub1(3:Hf, 2:Wf-1,:), 3) + ...
                              mean(sub1(2:Hf-1, 1:Wf-2,:) .* sub1(2:Hf-1, 3:Wf,:), 3)) / 2;
    cum2(2:Hf-1, 2:Wf-1) = (mean(sub2(1:Hf-2, 2:Wf-1,:) .* sub2(3:Hf, 2:Wf-1,:), 3) + ...
                              mean(sub2(2:Hf-1, 1:Wf-2,:) .* sub2(2:Hf-1, 3:Wf,:), 3)) / 2;
% cum1 = cum1.^0.85;
% cum2 = cum2.^0.85;
    % imwrite(uint16(cum1 ./ max(cum1(:)+eps) .* 65535), ...
    %     fullfile(out_dir, sprintf('SR%d_SOFI_ch1.tif', f)));
    % imwrite(uint16(cum2 ./ max(cum2(:)+eps) .* 65535), ...
    %     fullfile(out_dir, sprintf('SR%d_SOFI_ch2.tif', f)));

    %% 6. sparse_main 稀疏约束
    % sp1 = sparse_main(cum1.^0.5, sp_ch1_mu, sp_ch1_sigmat, sp_ch1_l1, sp_ch1_iter, sp_ch1_backg).^2;
    % sp2 = sparse_main(cum2.^0.5, sp_ch2_mu, sp_ch2_sigmat, sp_ch2_l1, sp_ch2_iter, sp_ch2_backg).^2;
sp1 = cum1;
sp2 = cum2;
    %% 7. 3× Fourier 上采样
    sp1_up = abs(fourierInterpolation(sp1.^0.5, [finter_sr, finter_sr], 'lateral')).^2;
    sp2_up = abs(fourierInterpolation(sp2.^0.5, [finter_sr, finter_sr], 'lateral')).^2;
    sp1_up = max(sp1_up, 0) ./ max(sp1_up(:) + eps);
    sp2_up = max(sp2_up, 0) ./ max(sp2_up(:) + eps);

    %% 8. 第二次 RL (PSF 按上采样倍数缩放)
    Ipsf_p1 = generate_psfv0(FWHM_post(1) * finter_sr);
    Ipsf_p2 = generate_psfv0(FWHM_post(2) * finter_sr);
    SACD1  = double(abs(deconvlucy(sp1_up, Ipsf_p1.^order(1), iter_post(1))));
    SACD2  = double(abs(deconvlucy(sp2_up, Ipsf_p2.^order(2), iter_post(2))));
    SACD1  = SACD1 ./ max(SACD1(:) + eps);  SACD1 = SACD1.^0.5;  SACD1 = SACD1 ./ max(SACD1(:) + eps);
    SACD2  = SACD2 ./ max(SACD2(:) + eps);  SACD2 = SACD2.^0.5;  SACD2 = SACD2 ./ max(SACD2(:) + eps);

    %% 累积多帧
    if f == 1
        accum_SACD1 = SACD1;
        accum_SACD2 = SACD2;
        accum_SOFI1 = sp1_up;
        accum_SOFI2 = sp2_up;
    else
        accum_SACD1(:,:,f) = SACD1;
        accum_SACD2(:,:,f) = SACD2;
        accum_SOFI1(:,:,f) = sp1_up;
        accum_SOFI2(:,:,f) = sp2_up;
    end

    % clear imgstack_f comp_f sig1 sig2 bg1 bg2 RL1 RL2 sub1 sub2 sp1 sp2 sp1_up sp2_up;
% end

%% ====== 裁剪 padding + 保存 ======
c1   = paddingfactor * finter_sr + 1;
crop = @(x) x(c1:end-paddingfactor*finter_sr, c1:end-paddingfactor*finter_sr, :);

imwritestack(uint16(crop(accum_SACD1) .* 65535), fullfile(out_dir, 'TMISACD_ch1.tif'));
imwritestack(uint16(crop(accum_SACD2) .* 65535), fullfile(out_dir, 'TMISACD_ch2.tif'));
imwritestack(uint16(crop(accum_SOFI1) .* 65535), fullfile(out_dir, 'SOFI_ch1.tif'));
imwritestack(uint16(crop(accum_SOFI2) .* 65535), fullfile(out_dir, 'SOFI_ch2.tif'));

fprintf('  [step3] SACD重建完成 → %s\n', out_dir);

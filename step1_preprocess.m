% =========================================================================
% step1_preprocess.m — 数据加载 + 预处理 (含文件循环)
% =========================================================================
% 每个文件的输出独立存放在: result/{exp_name}/{name}/
%   ├── {name}_raw.tif
%   ├── {name}_gaussRL.tif
%   ├── TMI_unmix/      (step2 输出)
%   ├── mask/           (step2b 输出)
%   └── SACD/           (step3 输出)
% =========================================================================

clc; close all;

addpath('TMI_script/TMI_script/TMISACDF');
addpath('TMI_script/TMI_script/sparse/F');

%% ====== 配置（从本机 tmisacd_config.m 加载，不提交）======
% 复制 tmisacd_config_template.m → tmisacd_config.m 并填入本机路径。
tmisacd_config;
% 注入: base_dir, file_pattern, name_pattern, idx_list,
%       exp_name, confname, rsFPs, frame, skip

result_root   = 'result';

%% ====== 预处理参数 ======
sigma         = [1, 1, 0];
FWHM_TMI      = 3;
iterTMI       = 3;
paddingfactor = 10;
RL1_on        = 0;

fprintf('处理 %d 个文件\n', numel(idx_list));

%% ====== 逐文件 ======
for idx = idx_list
    imgname = [base_dir, sprintf(file_pattern, idx)];
    fname   = sprintf(name_pattern, idx);
    fprintf('\n========== [%03d] %s ==========\n', idx, fname);

    % --- 每个文件独立的输出根目录 ---
    FILE_OUT_DIR = fullfile(result_root, exp_name, fname);
    if ~exist(FILE_OUT_DIR, 'dir'), mkdir(FILE_OUT_DIR); end
    fprintf('  输出目录: %s\n', FILE_OUT_DIR);

    % --- 加载 ---
    data_all = double(imreadstack(imgname));
    n_SR_frames = floor(size(data_all, 3) / skip);
    fprintf('  尺寸: %d x %d x %d, SR组数: %d\n', size(data_all, 1), size(data_all, 2), size(data_all, 3), n_SR_frames);
    imgstack = data_all(:, :, 1:frame);

    % --- 填充 + 归一化 ---
    imgstack = padarray(imgstack, [paddingfactor, paddingfactor, 0], 'symmetric');
    imgstack = imgstack ./ max(imgstack(:));
    imwritestack(uint16(imgstack .* 65535), fullfile(FILE_OUT_DIR, sprintf('%s_raw.tif', fname)));

    % --- 高斯去噪 ---
    if sigma(3) ~= 0
        imgstack_gauss = imgaussfilt3(imgstack, sigma);
    else
        imgstack_gauss = imgstack;
    end

    % --- RL 解卷积 ---
    IpsfTMI = generate_psfv0(FWHM_TMI);
    stack_RL = zeros(size(imgstack_gauss));
    for i = 1:frame
        stack_RL(:, :, i) = deconvlucy(imgstack_gauss(:, :, i), IpsfTMI, iterTMI);
    end

    if RL1_on == 1
        stack1_for_TMI = stack_RL;
    else
        stack1_for_TMI = imgstack_gauss;
    end

    imwritestack(uint16(stack_RL ./ max(stack_RL(:)) .* 65535), ...
        fullfile(FILE_OUT_DIR, sprintf('%s_gaussRL.tif', fname)));

    fprintf('  预处理完成\n');

    % ---- 链式运行后续步骤 ----
    run('step2_TMI_unmix');
    run('step2b_mask_postprocess');
    run('step3_SACD_reconstruct');
end

fprintf('\n===== 全部完成 =====\n');

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

addpath(genpath('F'));

%% ====== 配置（从本机 tmisacd_config.m 加载，不提交）======
% 复制 tmisacd_config_template.m → tmisacd_config.m 并填入本机路径。
tmisacd_config;
% 注入: base_dir, file_pattern, name_pattern, idx_list,
%       exp_name, confname, rsFPs, frame, skip

result_root   = 'result';

%% ====== 预处理参数 ======
FWHM_TMI      = 3.1;
iterTMI       = 5;
paddingfactor = 10;
RL1_on        = 1;

fprintf('处理 %d 个文件\n', numel(idx_list));

%% ====== 逐文件 ======
for idx = idx_list
    imgname = [base_dir, sprintf(file_pattern, idx)];
    fname   = sprintf(name_pattern, idx);

    % --- ROI（可选，由 roi_map 指定；无定义或为空则全图）---
    if exist('roi_map', 'var') && idx <= numel(roi_map) && ~isempty(roi_map{idx})
        roi   = roi_map{idx};
        fname = [fname, sprintf('_r%d-%d_c%d-%d', roi(1), roi(2), roi(3), roi(4))];
    else
        roi = [];
    end

    fprintf('\n========== [%03d] %s ==========\n', idx, fname);

    % --- 每个文件独立的输出根目录 ---
    FILE_OUT_DIR = fullfile(result_root, exp_name, fname);
    if ~exist(FILE_OUT_DIR, 'dir'), mkdir(FILE_OUT_DIR); end
    fprintf('  输出目录: %s\n', FILE_OUT_DIR);

    % --- 加载 ---
    data_all = double(imreadstack(imgname));
    if ~isempty(roi)
        data_all = data_all(roi(1):roi(2), roi(3):roi(4), :);
        fprintf('  ROI 裁剪: 行 %d-%d, 列 %d-%d  → %dx%d\n', ...
            roi(1), roi(2), roi(3), roi(4), roi(2)-roi(1)+1, roi(4)-roi(3)+1);
    end
    n_SR_frames = floor(size(data_all, 3) / skip);
    fprintf('  尺寸: %d x %d x %d, SR组数: %d\n', size(data_all, 1), size(data_all, 2), size(data_all, 3), n_SR_frames);

    % --- WF（全帧均值，无 padding）---
    wf = mean(data_all, 3);
    wf = wf ./ max(wf(:));
    imwrite(uint16(wf .* 65535), fullfile(FILE_OUT_DIR, 'WF.tif'));

    imgstack = data_all(:, :, 1:frame);

    % --- 填充 + 归一化 ---
    imgstack = padarray(imgstack, [paddingfactor, paddingfactor, 0], 'symmetric');
    imgstack = imgstack ./ max(imgstack(:));
    cp = 1 + paddingfactor;
    imwritestack(uint16(imgstack(cp:end-paddingfactor, cp:end-paddingfactor, :) .* 65535), ...
        fullfile(FILE_OUT_DIR, 'raw.tif'));

    % --- 解混前 RL 解卷积 (总开关 RL1_on) ---
    %   RL1_on=1: 解混前先做 RL → 再 PCA → 再解混; 这份 RL 数据同时传给 step3 时变分离。
    %   RL1_on=0: 不做 RL; step2 PCA / step3 时变分离都用未 RL 的数据。
    % stack1_for_TMI 即 "PCA 之前的时序数据"，step3 直接复用它做时变分离。
    if RL1_on == 1
        IpsfTMI = generate_psfv0(FWHM_TMI);
        stack_RL = zeros(size(imgstack));
        for i = 1:frame
            stack_RL(:, :, i) = deconvlucy(imgstack(:, :, i), IpsfTMI, iterTMI);
        end
        stack1_for_TMI = stack_RL;
        imwritestack(uint16(stack_RL ./ max(stack_RL(:)) .* 65535), ...
            fullfile(FILE_OUT_DIR, 'gaussRL.tif'));
    else
        stack1_for_TMI = imgstack;
    end

    fprintf('  预处理完成\n');

    % ---- 三种解混方法各跑完整 TMISACD 流程 (用 UNMIX_TAG 区分输出目录) ----
    % 输出: mask_unmix/SACD_unmix, mask_phasor/SACD_phasor, mask_rlsu/SACD_rlsu

    % 方法 1: fmincon sum-to-one 解混
    % run('step2_TMI_unmix');
    % run('step2b_mask_postprocess');
    % % run('step3_SACD_reconstruct');

    % 方法 2: Phasor 相量解混
    % run('step2_TMI_phasor');
    % run('step2b_mask_postprocess');
    % % run('step3_SACD_reconstruct');

    % 方法 3: RLSU (Richardson-Lucy) 解混
    run('step2_TMI_rlsu');
    run('step2b_mask_postprocess');
    run('step3_SACD_reconstruct');
end

fprintf('\n===== 全部完成 =====\n');

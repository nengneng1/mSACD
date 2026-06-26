function result = darksection_stack(data, params, background, pad, thres, divide, dark_blend, pad_size, dark3d)
% darksection_stack — Dark sectioning 背景去除 (共享函数, step2/step3 通用)
%   输入  data       : H×W×Nz, 归一化 [0,1]
%   输出  result     : 同尺寸, 去背景后并重新归一化到 [0,1]
%
% 算法: Kaiming He et al., "Single Image Haze Removal Using Dark Channel
%        Prior" (IEEE TPAMI 2011)。移植自 Caoruijie & Xipeng (PKU)。
%
% 参数:
%   params      : struct(NA, emwavelength, pixelsize, factor)
%   background  : 0=轻度(1次,deg=4,dep=2), 1=强力(2次,deg=[6,3,1.2],dep=[3,3,2])
%   pad         : 1=对称填充, 0=零填充
%   thres       : 背景阈值 (0–255 尺度, 典型 70)
%   divide      : Hi/Lo 分频点 (0.5)
%   dark_blend  : 去背景与原始的混合比 (1=完全去背景, 0=不变)
%   pad_size    : 边缘填充分母 (pad_amount = floor(Nx/pad_size)+1)
%   dark3d      : true =用时间均值估一张背景, 对所有帧统一减 (时间一致, 解混前必用)
%                 false=逐帧独立估计 (空间最优, 但会改变 decay 时序形状)
%
% 依赖: separateHiLo, confirm_block, dehaze_fast2 (及其子函数)
%       已复制到本 F 文件夹, 随脚本同目录调用即可。

% 缩放到 0-255
data = data * 255;

% 记录原始尺寸, 扩展为正方形 (separateHiLo 截止频率计算基于方形)
[Nx0, Ny0, Nz] = size(data);
if Ny0 > Nx0, data(Nx0+1:Ny0, :, :) = 0; end
if Ny0 < Nx0, data(:, Ny0+1:Nx0, :) = 0; end
[Nx, Ny, ~] = size(data);

% deg/dep/hl 矩阵
if background == 1
    maxtime    = 2;
    deg_matrix = [6, 3, 1.2];
    dep_matrix = [3, 3, 2];
    hl_matrix  = [1, 1, 1];
else
    maxtime    = 1;
    deg_matrix = [4];
    dep_matrix = [2];
    hl_matrix  = [1];
end

data_raw = data;
pad_r = floor(Nx / pad_size) + 1;
pad_c = floor(Ny / pad_size) + 1;
cx = floor(Nx/pad_size)+2 : floor(Nx/pad_size)+Nx+1;
cy = floor(Ny/pad_size)+2 : floor(Ny/pad_size)+Ny+1;

if dark3d
    % ---- 3D 模式: 用时序均值估背景, 对所有帧施加相同矫正 (保持 decay 形状) ----
    avg_data = mean(data, 3);
    for time = 1:maxtime
        deg = deg_matrix(time);  dep = dep_matrix(time);  hl = hl_matrix(maxtime);
        if pad == 1
            avg_padded = padarray(avg_data, [pad_r, pad_c], 'symmetric');
        else
            avg_padded = padarray(avg_data, [pad_r, pad_c]);
        end
        p = params;  p.Nx = size(avg_padded,1);  p.Ny = size(avg_padded,2);
        [Hi, Lo, lp, EL] = separateHiLo(avg_padded, p, deg, divide);
        block_size        = confirm_block(p, lp);
        Lo_proc           = dehaze_fast2(Lo, 0.95, block_size, EL, dep, thres);
        result_padded     = Lo_proc/hl + Hi;
        correction        = result_padded(cx, cy) - avg_data;   % 负值 = 减背景
        for jz = 1:Nz
            data(:,:,jz) = data(:,:,jz) + correction;
        end
        if time < maxtime, avg_data = mean(data, 3); end
    end
    result_final = data;

else
    % ---- 2D 模式: 逐帧独立估计 ----
    result_stack = zeros(Nx, Ny, Nz);
    if pad == 1
        img_pad_init = padarray(data(:,:,1), [pad_r, pad_c], 'symmetric');
    else
        img_pad_init = padarray(data(:,:,1), [pad_r, pad_c]);
    end
    image = zeros(size(img_pad_init,1), size(img_pad_init,2), Nz);
    for jz = 1:Nz
        if pad == 1
            image(:,:,jz) = padarray(data(:,:,jz), [pad_r, pad_c], 'symmetric');
        else
            image(:,:,jz) = padarray(data(:,:,jz), [pad_r, pad_c]);
        end
    end
    p = params;  p.Nx = size(image,1);  p.Ny = size(image,2);
    for time = 1:maxtime
        deg = deg_matrix(time);  dep = dep_matrix(time);  hl = hl_matrix(maxtime);
        for jz = 1:Nz
            [Hi, Lo, lp, EL] = separateHiLo(squeeze(image(:,:,jz)), p, deg, divide);
            block_size        = confirm_block(p, lp);
            Lo_proc           = dehaze_fast2(Lo, 0.95, block_size, EL, dep, thres);
            result_full       = Lo_proc/hl + Hi;          % 含 padding 的完整结果
            result_stack(:,:,jz) = result_full(cx, cy);   % 裁回原始尺寸
        end
        data = result_stack;
        for jz = 1:Nz
            if pad == 1
                image(:,:,jz) = padarray(data(:,:,jz), [pad_r, pad_c], 'symmetric');
            else
                image(:,:,jz) = padarray(data(:,:,jz), [pad_r, pad_c]);
            end
        end
    end
    result_final = result_stack;
end

% 混合 + 裁回原始非方形尺寸
result_final = dark_blend * result_final + (1 - dark_blend) * data_raw;
result_final = max(result_final, 0);
if Nx0 ~= Nx, result_final(Nx0+1:Nx, :, :) = []; end
if Ny0 ~= Ny, result_final(:, Ny0+1:Ny, :) = []; end

% 归一化回 [0,1]
result = result_final ./ max(result_final(:) + eps);
end

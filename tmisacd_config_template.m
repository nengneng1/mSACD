% TMISACD_CONFIG_TEMPLATE  TMI-SACD 双色流程本机路径与数据集配置模板
%
% 使用方法:
%   cp tmisacd_config_template.m tmisacd_config.m
%   然后修改 tmisacd_config.m 填入本机实际路径。
%   tmisacd_config.m 已加入 .gitignore，不会被提交。
%
% ══════════════════════════════════════════════════════════════════
% 填写指引
% ══════════════════════════════════════════════════════════════════
%
% 【base_dir】 tif 数据所在文件夹，末尾加 /
%   示例: 'D:/data/experiment/tif_even/'
%
% 【file_pattern】 文件名 sprintf 格式，参数: idx（整数序号）
%   示例: 'lenti_MAP4_MSG+Mito_rsfastlime_%03d_s1.tif'
%   验证: sprintf(file_pattern, 2) 应还原出第一个文件名
%
% 【name_pattern】 输出名格式，无空格无特殊字符，参数: idx
%   示例: 'lenti_MAP4_MSG_Mito_rsfastlime_%03d'
%
% 【idx_list】 要处理的文件序号，如 2:2:44（偶数批量）或 [2 4 6]（指定）
%
% 【confname】 荧光关闭曲线 .mat 文件名（不含扩展名），需放在工作目录下
% ══════════════════════════════════════════════════════════════════

%% ====== 路径（必填）======
base_dir     = '';   % tif 所在文件夹，末尾加 /
                     % 示例: 'D:/data/experiment/tif_even/'

file_pattern = '';   % 文件名 sprintf 格式，参数: idx
                     % 示例: 'lenti_MAP4_MSG+Mito_rsfastlime_%03d_s1.tif'

name_pattern = '';   % 输出名 sprintf 格式，参数: idx，无特殊字符
                     % 示例: 'lenti_MAP4_MSG_Mito_rsfastlime_%03d'

%% ====== 数据集配置 ======
idx_list = [];       % 文件序号，如 2:2:44（偶数批量）或 [2 4 6]（指定）

% ROI（按文件序号指定，不需要裁剪则留空 []）
% 格式: [r1, r2, c1, c2]（行起止, 列起止，MATLAB 1-indexed）
% 如何获取: FIJI 打开第1帧 → 矩形工具框选 → Image > Show Info 读取
%   BoundingBox: x, y, width, height
%   转换: [r1,r2,c1,c2] = [y+1, y+height, x+1, x+width]
% 输出名会自动附加 _r{r1}-{r2}_c{c1}-{c2} 后缀以区分不同裁剪区域。
roi_map = cell(1, max(idx_list(:)));
% roi_map{2}  = [r1, r2, c1, c2];   % 按需填写，其余自动全图

%% ====== 实验参数 ======
exp_name = '';       % 实验名，用作输出子目录名，示例: 'MAP4_MSG_Mito'
confname = '';       % 衰减曲线 .mat 文件名（不含扩展名），示例: 'tomm_map4'
rsFPs    = 2;        % 荧光团个数
frame    = 20;       % SACD 重建使用的帧数
skip     = 40;       % 每组超分辨帧之间的原始帧跳跃数

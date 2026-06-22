% 获取用户输入的图像文件夹路径
folder_path = input('请输入包含图像的文件夹路径: ', 's');

% 检查文件夹是否存在
if ~exist(folder_path, 'dir')
    error('指定的文件夹不存在');
end

% 获取用户输入的图像文件名
img_name1 = input('请输入第一张图像的文件名: ', 's');
img_name2 = input('请输入第二张图像的文件名: ', 's');
% C3_sn2n_1.tif
% 构建完整的图像路径
img_path1 = fullfile(folder_path, img_name1);
img_path2 = fullfile(folder_path, img_name2);

% 读取图像
try
    % 检查文件是否存在
    if ~exist(img_path1, 'file')
        error(['第一张图像 "', img_name1, '" 不存在']);
    end
    if ~exist(img_path2, 'file')
        error(['第二张图像 "', img_name2, '" 不存在']);
    end
    
    img1 = imread(img_path1);
    img2 = imread(img_path2);
    
    % 确保图像类型一致
    if ~isa(img2, class(img1))
        img2 = cast(img2, class(img1));
    end
    
    % 确保图像尺寸相同
    if size(img1) ~= size(img2)
        img2 = imresize(img2, size(img1));
        fprintf('警告: 两张图像尺寸不同，已将第二张图像调整为与第一张相同的尺寸\n');
    end
%     C3_nor-sn2n.tif
    % 使用MATLAB自带的psnr函数计算
    psnr_value = psnr(img1, img2);
    ssim_value = ssim(img1, img2);
    
    % 显示结果
    fprintf('图像 "%s" 和 "%s" 的PSNR值为: %.4f dB\n', ...
        img_name1, img_name2, psnr_value);
    fprintf('图像 "%s" 和 "%s" 的SSIM值为: %.4f \n', ...
        img_name1, img_name2, ssim_value);
catch ME
    fprintf('发生错误: %s\n', ME.message);
end    
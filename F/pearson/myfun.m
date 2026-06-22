function myfun(refImagePath, procFolder)
    % 读取固定的参考图像
    refImage = imread(refImagePath);
    
    % 确保参考图像是灰度图像
    if ndims(refImage) > 2
        error('The reference image must be a grayscale image.');
    end
    
    % 获取处理后图像列表
    procImages = dir(fullfile(procFolder, '*.tif'));  % 假设图像格式为.tiff
    
    % 初始化存储PSNR和SSIM值的数组
    numImages = length(procImages);
    psnrValues = zeros(numImages, 1);
    ssimValues = zeros(numImages, 1);
    
    % 遍历所有处理后的图像并计算PSNR和SSIM
    for i = 1:numImages
        % 读取处理后的图像
        procImagePath = fullfile(procFolder, procImages(i).name);
        procImage = imread(procImagePath);
        
        % 确保处理后的图像是灰度图像
        if ndims(procImage) > 2
            warning(['Image is not grayscale: ' procImagePath]);
            continue;
        end
        
        % 确保两幅图像具有相同的尺寸
        if ~isequal(size(refImage), size(procImage))
            warning(['Image sizes do not match: ' refImagePath ' and ' procImagePath]);
            continue;
        end
        
        % 计算PSNR
        psnrValues(i) = psnr(procImage, refImage);
        
        % 计算SSIM
        [ssimValues(i), ~] = ssim(im2double(procImage), im2double(refImage));
        
        % 显示进度
        disp(['Processed image: ', procImages(i).name]);
    end

    % 绘制PSNR和SSIM曲线
    figure;
    subplot(2, 1, 1);
    plot(1:numImages, psnrValues, '-o', 'LineWidth', 2);
    title('PSNR Values');
    xlabel('Image Index');
    ylabel('PSNR (dB)');
    grid on;
    
    subplot(2, 1, 2);
    plot(1:numImages, ssimValues, '-o', 'LineWidth', 2);
    title('SSIM Values');
    xlabel('Image Index');
    ylim([0 1]);
    ylabel('SSIM');
    grid on;
    
    % 保存结果到CSV文件
%     results = table((1:numImages)', psnrValues, ssimValues, ...
%                     'VariableNames', {'ImageIndex', 'PSNR', 'SSIM'});
%     writetable(results, 'psnr_ssim_results.csv');
end

% 调用函数，指定固定的参考图像路径和处理后图像的文件夹路径
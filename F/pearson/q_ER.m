clc
clear
folderPath = 'e';

tifFiles_red = dir(fullfile(folderPath, '*C561.tif'));
fileNames_red = {tifFiles_red.name};
fileNames_red = string(fileNames_red');
%% 
pearsonVal_all = zeros(1,37); 
ssimVal_all = zeros(1,37);
psnrVal_all = zeros(1,37);
pearsonVal_all2 = zeros(1,37);
ssimVal_all2 = zeros(1,37);
psnrVal_all2 = zeros(1,37);
figure('Position', [100, 100, 1200, 500]);  % 设置大尺寸画布
% stack = [9,13,14,15,18,19,20,22,24]
stack = [25]
%% mSACD
for n = stack
    imgname = fileNames_red(n);
    redRAW = imreadstack([imgname]);
    imTruth = redRAW(:,:,1);
    mSACD = imreadstack(['result/20250724e/',num2str(n),'_frame20/','RLon1_sigma110_fwhm2.52.5_iter3_3_finter4gaus11_order1_0.7_TMISACD_ch1.tif']);
    imSplit = mSACD(:,:,1);
    imSplit= imresize(imSplit,[0.25*size(mSACD,1),0.25*size(mSACD,2)],'cubic');
    imSplit = imSplit./max(imSplit(:));
    imSplit = imgaussfilt(imSplit,2);
    imSplit = imSplit./max(imSplit(:));

    imTruth = imTruth./max(imTruth(:));
    [optimizer, metric] = imregconfig('monomodal'); % 配置单模态配准参数
    imSplit= imregister(imSplit,imTruth,'translation', optimizer, metric, ...
        'DisplayOptimization', false);
    [pearsonVal, ssimVal, psnrVal] = calculateMetrics(imSplit, imTruth)
    pearsonVal_all(:,n) = pearsonVal;
    ssimVal_all(:,n) = ssimVal;
    psnrVal_all(:,n) = psnrVal;
end
%% TMI
for n =stack
    imgname = fileNames_red(n);
    redRAW = imreadstack([imgname]);
    imTruth = redRAW(:,:,1);
    TMI = imreadstack(['result/20250724e/',num2str(n),'_frame20/','TMIch1.tif']);
    imSplit2 = TMI(:,:,1);
    imSplit2 = imSplit2./max(imSplit2(:));
    % imSplit2 = imgaussfilt(imSplit2,2);
    imTruth = imTruth./max(imTruth(:));
    [optimizer, metric] = imregconfig('monomodal'); % 配置单模态配准参数
    imSplit2= imregister(imSplit2,imTruth,'translation', optimizer, metric, ...
        'DisplayOptimization', false);
    [pearsonVal2, ssimVal2, psnrVal2] = calculateMetrics(imSplit2, imTruth)
    pearsonVal_all2(:,n) = pearsonVal2;
    ssimVal_all2(:,n) = ssimVal2;
    psnrVal_all2(:,n) = psnrVal2;
end
%% 
subplot(121);
% 合并两种算法的Pearson数据
pearsonData = [pearsonVal_all; pearsonVal_all2]'; 
boxplot(pearsonData, 'Labels', {'mSACD', 'TMI'});
title('Pearson系数对比');
ylabel('相关系数');
grid on;

subplot(122);
% 合并两种算法的SSIM数据
ssimData = [ssimVal_all; ssimVal_all2]'; 
boxplot(ssimData, 'Labels', {'mSACD', 'TMI'});
title('SSIM相似度对比');
ylabel('结构相似度');
grid on;

%% 
% 生成时间戳文件名
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('Algorithm_Comparison_%s.png', timestamp);

% 设置保存参数
set(gcf, 'PaperPositionMode', 'auto'); 
print(gcf, filename, '-dpng', '-r600');  % 600dpi高清输出
fprintf('对比图已保存为: %s\n', filename);

%% 统计值表格

statsTable = table(...
    mean(pearsonVal_all'), std(pearsonVal_all'), ...
    mean(ssimVal_all'), std(ssimVal_all'), ...
    mean(pearsonVal_all2'), std(pearsonVal_all2'), ...
    mean(ssimVal_all2'), std(ssimVal_all2'), ...
    'VariableNames', ...
    {'TMISACD_Pearson_Mean', 'mSACD_Pearson_STD', ...
     'TMISACD_SSIM_Mean', 'mSACD_SSIM_STD', ...
     'TMI_Pearson_Mean', 'TMI_Pearson_STD', ...
     'TMI_SSIM_Mean', 'TMI_SSIM_STD'});

disp('算法性能统计:');
disp(statsTable);
writetable(statsTable, 'Performance_Statistics.xlsx');
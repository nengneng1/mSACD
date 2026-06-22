clc;
clear;
mean_values = zeros(10,1);
for i = 1:7
img_TMI = double(imreadstack(['TMI',num2str(i),'.tif']));
img_gt = double(imreadstack(['GT',num2str(i),'.tif']));

img_TMI = img_gt (:,:,2);

img_gt = img_gt (:,:,1);
img_gt = imresize(img_gt,[size(img_TMI,1),size(img_TMI,2)]);
img_TMI = img_TMI./max(max(img_TMI(:)));
img_gt = img_gt./max(max(img_gt(:)));
    % 计算每个图像的平均值
    mean_imgTMI = mean(img_TMI(:));
    mean_img_gt = mean(img_gt(:));

    % 计算协方差和标准差
    covar = mean((img_TMI(:) - mean_imgTMI) .* (img_gt(:) - mean_img_gt));
    std_dev_img1 = sqrt(mean((img_TMI(:) - mean_imgTMI).^2));
    std_dev_img2 = sqrt(mean((img_gt(:) - mean_img_gt).^2));

    % 避免除以零的情况
    if std_dev_img1 == 0 || std_dev_img2 == 0
        error('一个或两个图像的标准偏差为零，无法计算皮尔森相关系数');
    end

    % 计算皮尔森相关系数
    r = covar / (std_dev_img1 * std_dev_img2);
data(i,1) = r;
end

mean_values = mean(data);

std_devs = std(data);

% 计算每个条件的标准误差
n = size(data, 1); % 假设每列的数据量相同
std_errors = std_devs / sqrt(n);
% 
% % 示例数据
% mean_values = [0.85, 0.82, 0.88, 0.86, 0.84, 0.83, 0.87, 0.85, 0.86];
% std_errors = [0.02, 0.03, 0.01, 0.02, 0.03, 0.02, 0.01, 0.02, 0.03];
% 
% 条形图的颜色
colors = {'m', 'w', 'k', 'm', 'w', 'k', 'g', 'w', 'k'};

% 绘制条形图
figure;

b =bar(1,mean_values,0.01);

 set(b(1),'facecolor',[0,0.5,0.5])
hold on;

% 添加误差线
for i = 1:length(mean_values)
    errorbar(i, mean_values(i), std_errors(i), 'LineStyle', 'none', 'Marker', '.', 'Color', colors{i});
end

% 设置图表属性
% xlabel('Conditions');
ylabel('Pearson''s correlation');
title('Pearson''s correlation for different channel');
% legend('channel 1', 'Condition 2', 'Condition 3', 'Condition 4', 'Condition 5', 'Condition 6', 'Condition 7', 'Condition 8', 'Condition 9');
xticklabels({'channel 1', 'YFP-FLAG with the other five non-Flag constructs', 'Dronpa-FLAG only', 'Dronpa-FLAG with the other five non-Flag constructs', 'Skyylan62A-FLAG only', 'Skyylan62A-FLAG with the other five non-Flag constructs', 'rsFastlime-FLAG only', 'rsEGFP2-E-FLAG only', 'rsGreenF-E-FLAG only'});
% xtickangle(45);
% 添加数据点
for i = 1:size(data, 1)
    for j = 1:size(data, 2)
        scatter(1, data(i, j), 10, 'filled', 'Color', colors{i});
    end
end
% 添加显著性标记
% text(1.5, 0.9, 'n.s.', 'HorizontalAlignment', 'center');
% text(4.5, 0.9, 'n.s.', 'HorizontalAlignment', 'center');
% text(7.5, 0.9, 'n.s.', 'HorizontalAlignment', 'center');
% text(10.5, 0.9, 'n.s.', 'HorizontalAlignment', 'center');
% text(13.5, 0.9, 'n.s.', 'HorizontalAlignment', 'center');

hold off;
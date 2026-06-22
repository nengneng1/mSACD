clc;clear;close all;
addpath('F');
% initialization

folder = '202504122_3colors';

mkdir('result');
num = 85;
channel = [405,561,640];%%
channel = channel(3);
imgname = ['tomm20-sky62A&Mito-rstagRFP&KDEL-Dsred-60X_20%_488_40%_561__',num2str(num),'.vsi - C',num2str(channel)];

imgfolder = ['result/',folder,'/',num2str(num),'sparse'];mkdir(imgfolder);

%% parmter
finter = 2;%上采样倍数

fidelity0 = 50 ;%保真度，越低则越平滑
fidelity_z0 = 0; %z轴保真度，各项同性默认为1
sparsity0 = 3; %稀疏度，提高分辨率、去除离焦信号
backg0 = 3;%背景参数：当设置为0时，不滤除背景；当设置为大于1的数时，越接近1滤除效果越明显

iter = 5;%解卷积迭代次数
pixel = 1*10^-9;%像素尺寸
FWHM = 3.5* 10^-9;%卷积核半高宽

%% load
Ipsf2D = PSFget2D((finter*FWHM/pixel));
data = imreadstack([imgname,'.tif']);
data  = data (:,:,1);
imWS(data,[imgfolder,'/raw']);
%% main
data_upsamp = abs(single(fourierInterpolation(data,[finter,finter],'lateral')));%三维数据加1维的finter
data_upsamp = data_upsamp./max(max(data_upsamp(:)));
% imWS(data_upsamp,[imgfolder,'/fourier']);

[data_sparse] = sparse_main(abs(data_upsamp),fidelity0,fidelity_z0,sparsity0,100,backg0);
data_sparse = data_sparse./max(max(data_sparse(:)));
% imWS(data_sparse,[imgfolder,'/','sparse']);

data_RL = RL3D(data_sparse,Ipsf2D,iter,1);
data_RL = abs(data_RL);
% imWS(data_RL,[imgfolder,'/','deconv']);


data_all(:,:,1) = data_upsamp;
data_all(:,:,2) =data_sparse;
data_all(:,:,3) = data_RL;
imWS(data_all,[imgfolder,'/','all']);
function output = percennorm(data, miper, maper)
if nargin < 3, maper = 100; end
if nargin < 2, miper = 0; end
% data = single(data);
max(data(:));
datamin = prctile(data(:), miper);
datamax = prctile(data(:), maper);
output = (data - datamin) / (datamax - datamin);
output(output > 1) = 1;
output(output < 0) = 0;
end

function [output] = weightmat(input)
input(input<0.1)=0.1;
output = ones(size(input))./input;
output = output./mean(mean(mean(output)));
end

function  [output] = PSFget2D(sigma)
%sigma_3 = sigma;
sigma_2 = [sigma sigma];
psfsigma = sigma_2;
psfN = ceil(psfsigma./sqrt(8*log(2)) * sqrt(-2 * log(0.0002))) + 1;
psfN = psfN * 2 + 1;
output = Gauss(psfsigma,psfN);
end
function [data_decon] = deconv_main2(fornow,pixel_rescale,sigma,iter,gpuswitch)
%% parameter selection
sigma_2 = [sigma sigma sigma];
%% PSF
psfsigma = sigma_2./pixel_rescale;
psfN = ceil(psfsigma./sqrt(8*log(2)) * sqrt(-2 * log(0.0002))) + 1;
psfN = psfN * 2 + 1;
psf = Gauss(psfsigma,psfN);
%% LR deconv
fornow2(:,:,3:size(fornow,3)+2) = fornow(:,:,:);
fornow2(:,:,1) = fornow(:,:,1);
fornow2(:,:,2) = fornow(:,:,1);
data_decon0 = RL3D(fornow2,psf.^4,iter,gpuswitch);
data_decon = data_decon0(:,:,3:size(data_decon0,3));
%% save data
data_decon = single(gather(data_decon));
data_decon = data_decon./max(max(max(data_decon)));
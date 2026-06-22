clear all
close all
clc
eps=0.000001;
filt_no_PC=3; %3
%%
temp=double(imreadstack('ERtomm.tif'));
temp_bin_corr=temp;
NoiseCorrection=sqrt(mean(temp_bin_corr,[1 2])+eps);
temp_norm = temp_bin_corr ./ NoiseCorrection;
temp_norm(temp_norm<0)=0;
%%
vector_conv=((reshape(temp_norm,[size(temp_norm,1)*size(temp_norm,2) size(temp_norm,3)])));
%% calculate the data using PCA
[flimcoeff,scores,flimcoeffVar] = pca(vector_conv);
%%
%%reconstruct using PCA
filt_vector_PCA=(scores(:,1:filt_no_PC)*flimcoeff(:,1:filt_no_PC)')+mean(vector_conv,1);
filt_img_PCA=reshape(filt_vector_PCA,[size(temp_norm,1) size(temp_norm,2) size(temp_norm,3)]);
filt_img_PCA_conv=(filt_img_PCA.*NoiseCorrection(:,:,:));
filt_img_PCA_single=single(filt_img_PCA_conv);
imwritestack(filt_img_PCA_single./max(filt_img_PCA_single(:)).*65535,'denoised_3.tif')
imwritestack(temp./max(temp(:)).*65535,'raw.tif')

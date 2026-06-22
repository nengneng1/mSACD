clear;
clc;
frame = 30;
n = 103;
finter = 1;

folder = '20250321_TOMM';
imgfolder = ['result/',folder,'/',num2str(n),'_frame',num2str(frame),''];mkdir(imgfolder);
load('103_RLon0_sigma0.50.530_mask')
stack4SACD = imreadstack('103_RLon0raw_RL.tif');
stack4TMI = imreadstack('103_RLon0gauss_RL.tif');

composition_image= composition_image./max(composition_image(:));
imwritestack(double(composition_image .*65535),[imgfolder,'/','_mask.tif']);

composition_image_sparse = zeros(size(composition_image,1),size(composition_image,2),size(composition_image,3));
for i = 1:2
    composition_image_sparse(:,:,i) = sparse_main(composition_image(:,:,i),1000,0,1,100,5);
end
composition_image_sparse(:,:,3) = composition_image(:,:,3);
imwritestack(double(composition_image_sparse .*65535),[imgfolder,'/','_sparsemask.tif']);

composition_image_resize = zeros(finter*size(composition_image,1),finter*size(composition_image,2),size(composition_image,3));
composition_image_resize = abs(fourierInterpolation(composition_image_sparse ,[finter,finter,1],'lateral'));

% composition_imagel = composition_image_resize;

% middle1 = composition_imagel(:,:,1);
% middle2 = composition_imagel(:,:,2);
% middle3 = composition_imagel(:,:,3);
% middle2(middle1 > middle2) = middle2(middle1 > middle2)/par1;
% middle2(middle1 > middle2) = 0;
% middle2(middle2 > middle1) = 1;
% middle1(middle2 > middle1) = middle1(middle2 > middle1)/par1;
% middle1(middle2 > middle1) = 0;
% middle1(middle1 > middle2) = 1;
% composition_image(:,:,1) =  middle1;
% composition_image(:,:,2) =  middle2;


order = 2;
stack4SACD = percennorm(stack4SACD, 0,100);
stack4TMI = percennorm(stack4TMI, 0,100);
stack4SACD_resize = zeros(finter*size(stack4SACD,1),finter*size(stack4SACD,2),size(stack4SACD,3));
stack4SACD_resize = abs(fourierInterpolation(stack4SACD,[finter,finter,1],'lateral')); % fourier upsampling
stack4TMI_resize = zeros(finter*size(stack4TMI,1),finter*size(stack4TMI,2),size(stack4TMI,3));
stack4TMI_resize = abs(fourierInterpolation(stack4TMI,[finter,finter,1],'lateral')); % fourier upsampling
for ga = [0.5]
    for mul =[0.5];
        lag = 1;
        % for o = [0.3,0.4,0.5,0.6,0.7,1]
        FWHM2 = 3;
        iter2 = 5;%iteration times of post decovnolution
        Ipsf2 = generate_psfv0(FWHM2*finter);
        Ipsf_mito = generate_psfv0(2.5*finter);
        sparse_on = 1;
        gamma2 = ga;
        normax = 100;
        %% 2 sparse param
        fidelity0 = 200;%保真度，越低则越平滑
        fidelity_z0 = 0; %z轴保真度，各项同性默认为1
        sparsity0 = 0.1; %稀疏度
        back_sparse =0;%背景参数：当设置为0时，不滤除背景；当设置为大于1的数时，越接近1滤除效果越明显
        %% SACD cumulant

        maskimg1 = (composition_image_resize(:,:,1).*stack4SACD_resize.^(gamma2)).^(1/gamma2);
        maskimg2 = (composition_image_resize(:,:,2).*stack4SACD_resize.^(gamma2)).^(1/gamma2);
        maskimg3 = (composition_image_resize(:,:,3).*stack4SACD_resize.^(gamma2)).^(1/gamma2);
        maskimg1 = maskimg1./max(max(maskimg1(:)));
        maskimg2 = maskimg2./max(max(maskimg2(:)));
        maskimg3 = maskimg3./max(max(maskimg3(:)));
        % maskimgTMI1 = (composition_image(:,:,1).*stack4TMI_resize.^(gamma2)).^(1/gamma2);
        % maskimgTMI2 = (composition_image(:,:,2).*stack4TMI_resize.^(gamma2)).^(1/gamma2);
        % maskimgTMI3 = (composition_image(:,:,3).*stack4TMI_resize.^(gamma2)).^(1/gamma2);
        % maskimgTMI1 = maskimgTMI1./max(max(maskimgTMI1(:)));
        % maskimgTMI2 = maskimgTMI2./max(max(maskimgTMI2(:)));
        % maskimgTMI3 = maskimgTMI3./max(max(maskimgTMI3(:)));
        imwritestack(double(maskimg1.*65535),[imgfolder,'/','_mask1.tif']);
        imwritestack(double(maskimg2.*65535),[imgfolder,'/','_mask2.tif']);
        imwritestack(double(maskimg3.*65535),[imgfolder,'/','_mask3.tif']);
        %     imwritestack(double(maskimgTMI1.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'_par',num2str(par1*10),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2),num2str(FWHM4gauss3)],'_mul',num2str(mul),'_lag',num2str(lag),'_dif',num2str(diff_param),'_order',num2str(gamma_1),'_',num2str(gamma2),'_frame',num2str(frame),'_maskTMI1.tif']);
        % imwritestack(double(maskimgTMI2.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'_par',num2str(par1*10),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2),num2str(FWHM4gauss3)],'_mul',num2str(mul),'_lag',num2str(lag),'_dif',num2str(diff_param),'_order',num2str(gamma_1),'_',num2str(gamma2),'_frame',num2str(frame),'_maskTMI2.tif']);
        % imwritestack(double(maskimgTMI3.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'_par',num2str(par1*10),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2),num2str(FWHM4gauss3)],'_mul',num2str(mul),'_lag',num2str(lag),'_dif',num2str(diff_param),'_order',num2str(gamma_1),'_',num2str(gamma2),'_frame',num2str(frame),'_maskTMI3.tif']);
        stacksub1 = abs(maskimg1 - mul * mean(maskimg1,3));%获得涨落信号
        stacksub2 = abs(maskimg2 - mul * mean(maskimg2,3));%获得涨落信号
        stacksub3 = abs(maskimg3 - mul * mean(maskimg3,3));%获得涨落信号
        stacksub1 =  stacksub1(:,:,1:10);
        stacksub2 =  stacksub2(:,:,1:20);
        stacksub3 =  stacksub3(:,:,1:30);
        Nx=size(stacksub1,1);
        Ny=size(stacksub1,2);
        cumimg1=zeros(Nx,Ny);
        cumimg2=zeros(Nx,Ny);
        cumimg3=zeros(Nx,Ny);
        cumimg1(2:Nx-1,2:Ny-1) = (mean(stacksub1(1:Nx-2,2:Ny-1,:).*stacksub1(3:Nx,2:Ny-1,:),3) ...
            + mean(stacksub1(2:Nx-1,1:Ny-2,:).*stacksub1(2:Nx-1,3:Ny,:),3))./2;
        cumimg2(2:Nx-1,2:Ny-1) = (mean(stacksub2(1:Nx-2,2:Ny-1,:).*stacksub2(3:Nx,2:Ny-1,:),3) ...
            + mean(stacksub2(2:Nx-1,1:Ny-2,:).*stacksub2(2:Nx-1,3:Ny,:),3))./2;
        cumimg3(2:Nx-1,2:Ny-1) = (mean(stacksub3(1:Nx-2,2:Ny-1,:).*stacksub3(3:Nx,2:Ny-1,:),3) ...
            + mean(stacksub3(2:Nx-1,1:Ny-2,:).*stacksub3(2:Nx-1,3:Ny,:),3))./2;
        %     cumimg1  = cumulant0(stacksub1,2,lag);
        % cumimg2  = cumulant0(stacksub2,2,lag);


        clear maskimg3
        clear maskimg2
        clear maskimg1
        clear maskimgTMI1
        clear maskimgTMI2
        clear maskimgTMI3
        clear stacksub3
        clear stacksub2
        clear stacksub1
        %% 每个通道单独调节解卷积参数
        %1
        sparse1 = sparse_main(cumimg1.^(1/2),100,0,5,100,back_sparse);
        SACD1 = double(abs(deconvlucy(sparse1.^2, Ipsf2.^order, 3))); %the second RL deconvolution
        SACD1 = double(SACD1./max(SACD1(:)));
        SACD1 = SACD1.^(1/2);
        %2
        sparse2 = sparse_main(cumimg2.^(1/2),1000,0,2,100,back_sparse);
        SACD2 = double(abs(deconvlucy(sparse2.^2, Ipsf_mito.^order, 5))); %the second RL deconvolution
        SACD2 = double(SACD2./max(SACD2(:)));
        SACD2 = SACD2.^(1/2);
        %3
        sparse3 = sparse_main(cumimg3.^(1/2),500,0,2,100,0);
        SACD3 = double(abs(deconvlucy(sparse3.^2, Ipsf2.^order, 10))); %the second RL deconvolution
        SACD3 = double(SACD3./max(SACD3(:)));
        SACD3 = SACD3.^(1/2);

        %%
        [x,y] = size(SACD1)
        SACD2 = imresize(SACD2,[x,y]);
        SACD3 = imresize(SACD3,[x,y]);
        SACD_composi(:,:,1) =  SACD1;
        SACD_composi(:,:,2) =  SACD2;
        SACD_composi(:,:,3) =  SACD3;
        imwritestack(double(SACD_composi.*65535),[imgfolder,'/','_TMISACD.tif']);
    end

end

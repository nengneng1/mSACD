%%
%***************************************************************************
% TMI-SACD demo version 0.2
%***************************************************************************
%% initiation
clc;clear;close all;
folder = '20250506_3colors';
gamma_on = 1;
for n  = [1:2:37]
    for gamma_1 =[1]
        RL1_on = 0;
        sparseTMIon = 1;
        filename_front = 'ER-mstaygold&lenti-tomm20-sky62A_D2_60%_20ms_';
        num = num2str(n);
        filename_last = '.vsi - C488.tif';
        mkdir('result');
        imgname = [filename_front,num,filename_last];
        finter = 2;
        confname = '20ms_ER_tomm';%荧光关闭曲线文件名
        rsFPs = 2;%荧光团个数
        frame = 20; %frame number used for SACD reconstruction
        FWHM = 3.5;
        FWHM_TMI = 3.5;
        iter1 = 5; %iteration times of pre decovnolution
        iterTMI = 10;

        ifsub = 0;%是否减背景
        backg = 3;%减背景系数，越大减的越少
        baseline = 0;%图像基线信号，直接在图像数据上减除
        %% 
        sigma = [1,1,0];%对raw高斯去噪
        sigma2 = 0;%对衰减信号高斯去噪
        %% 
        par1 = 1;
        TMIon = 1;%是否进行通道拆解，如果设为0就只输出SACD结果
imgfolder = ['result/',folder,'/',num,'_frame',num2str(frame),''];mkdir(imgfolder);
        %% SACD RL1
        Ipsf1 = generate_psfv0(FWHM);
        IpsfTMI = generate_psfv0(FWHM_TMI);

        imgstack = double(imreadstack([imgname]));
        imgstack = imgstack(:,:,1:frame);
        imgstack = imgstack - baseline;
        imgstack(imgstack < 0) = 0;
        imwritestack(double(abs(imgstack)),[imgfolder,'/',num,'_raw.tif']);
        BACK = zeros([size(imgstack,1),size(imgstack,2),size(imgstack,3)]);
        if backg ~= 0%减背景
            for i =  1:frame
                BACK(:,:,i) = background_estimation(imgstack(:,:,i)./backg,1,6,'db6',3);
                BACK(BACK < 0) = 0;
            end
        else
            BACK = zeros([size(imgstack,1),size(imgstack,2)]);
        end
        if ifsub == 1
            imgstack = imgstack - BACK;
            imgstack(imgstack < 0) = 0;
        end

        imwritestack(double(abs(imgstack)),[imgfolder,'/',num,'_back_remove.tif']);
        % imwritestack(double(abs(BACK)),[imgfolder,'/',num,'_back1.tif']);
        % imwritestack(double(abs(BACK2)),[imgfolder,'/','_back2.tif']);
        stack4SACD = zeros([size(imgstack,1),size(imgstack,2),frame]);
        imgstack_gauss_RL1 = zeros([size(imgstack,1),size(imgstack,2),frame]);
        if sigma(3) ~=0
            imgstack_gauss =imgaussfilt3(imgstack,sigma);
        end
        if sigma(3) == 0
            imgstack_gauss = imgstack;
        end

        for i = 1:frame
            % stack(:,:,i) = sparse_main(imgstack(:,:,i),500,0,10,100,0); % the first RL deconvolution
            imgstack_gauss_RL1(:,:,i) = deconvlucy(imgstack_gauss(:,:,i),IpsfTMI,iterTMI); % the first RL deconvolution
        end
        if RL1_on == 1
            stack1 = imgstack_gauss_RL1;
        else
            stack1 = imgstack_gauss;
        end
        clear imgstack_gauss;

        if sparseTMIon ==1
            stack1 = sparse_main(stack1,200,0.5,0.1,100,0); % the first RL deconvolution
            stack1 = deconvlucy(stack1,IpsfTMI,iterTMI);
            imwritestack(stack1./max(stack1(:)).*65535,[imgfolder,'/',num,'_sparsestack.tif']);
        end

        if TMIon == 1
            %% TMI
            load([confname,'.mat']);

            if gamma_on == 1
                stack1 = stack1.^gamma_1;

                ch1 = ch1.^gamma_1;
                ch1 =ch1./max(ch1);
                ch2 = ch2.^gamma_1;
                ch2 =ch2./max(ch2);

            end
            plot(ch1);hold on
            plot(ch2);

            savefig([imgfolder,'/',num,'gamma',num2str(gamma_1),'_curve.fig'])
            % close all;
            imageWidth = size(stack1,1);
            imageHeight = size(stack1,2);
            imagedepth =  size(stack1,3);
            composition_imageraw1 = zeros(imageWidth, imageHeight, rsFPs);
            status_indicator1 = zeros(imageWidth, imageHeight) - 10;
            options = optimoptions('fmincon', 'Display', 'none','MaxIter', 5000);
            % options = optimoptions('fmincon', 'Display', 'none');
            stack1 = double(stack1./max(stack1(:)));
            % imwritestack(double(abs(stack1.*65535)),[imgfolder,'/',num,'_stack1.tif']);
            parfor ii = 1 : imageWidth
                for jj = 1: imageHeight
                    % if stack1(ii,jj,1) < (1 - ifsub) * double(BACK1(ii,jj)) || stack1(ii,jj,1) < 10/255 %
                    if stack1(ii,jj,1) < 500/65535 %
                        % ignore background pixels for faster computation
                        composition_imageraw1(ii,jj,:) = 0;
                    else
                        curve4s = squeeze(double(stack1(ii,jj,:)));
                        if sigma2 ~= 0
                            curve4s = imgaussfilt(curve4s, sigma2);
                        end
                        [predicted_vals1, ~, exitflag1, ~] = fmincon(@(coeffs)signalseparation2fmincon(ch1(1:frame), ch2(1:frame), curve4s, coeffs), ... % obj. fun.
                            [0.5,0.5], ... % initial guess
                            [], [], ... % inequalities (N/A)
                            ones(1,rsFPs), ... % equality constraint left hand side
                            1, ... % equality constraint right hand side -- this ensures they add up to 1
                            zeros(1,rsFPs), ... % lower bounds
                            ones(1,rsFPs), ... % upper bounds
                            [], ... % nonlinear options (N/A)
                            options); % see options variable above --  currently only disables print output
                        composition_imageraw1(ii,jj,:) = predicted_vals1(1:rsFPs);
                        status_indicator1(ii, jj)  = exitflag1; % we need this value to know if the optimization succeeded
                    end
                    disp(['Done with ' num2str(ii) ' of ' num2str(imageHeight)]);
                end
            end
            clear stack1;

            for i = 1:frame
                stack4SACD(:,:,i) = deconvlucy(imgstack(:,:,i),Ipsf1,iter1); % the first RL deconvolution
            end
            stack4SACD = percennorm(stack4SACD, 0,100);
            stack4SACD_resize = zeros(finter*size(stack4SACD,1),finter*size(stack4SACD,2),size(stack4SACD,3));
            stack4SACD_resize = abs(fourierInterpolation(stack4SACD,[finter,finter,1],'lateral')); % fourier upsampling

            imwritestack(double(abs(stack4SACD.*65535)),[imgfolder,'/',num2str(n),'_RLon',num2str(RL1_on),'raw_RL.tif']);
            imwritestack(double(abs(imgstack_gauss_RL1./max(imgstack_gauss_RL1(:)).*65535)),[imgfolder,'/',num2str(n),'_RLon',num2str(RL1_on),'gauss_RL.tif']);

            % for g  = [2]
            FWHM4gauss1 = 2;%用于模糊TMImask的高斯半高宽，因为拆解出的图像比较断裂，需要做高斯模糊来弥补一下
            FWHM4gauss2 = 2;


            %%
            composition_image = zeros(imageWidth, imageHeight, rsFPs);
            composition_image = min(max(composition_imageraw1, 0),1);
            % composition_image (composition_image < 0.4) = 0;
            %mask_gauss
            if FWHM4gauss1~=0
                composition_image(:,:,1) = imgaussfilt(composition_image(:,:,1), FWHM4gauss1);
            end
            if FWHM4gauss2~=0
                composition_image(:,:,2) = imgaussfilt(composition_image(:,:,2), FWHM4gauss2);
            end
            mask_mat_name = [imgfolder,'/',num2str(n),'_RLon',num2str(RL1_on),'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_mask.mat']
            save(mask_mat_name,'composition_image');
            imwritestack(double(composition_image .*65535),[imgfolder,'/',num2str(n),'_RLon',num2str(RL1_on),'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_mask.tif']);
            % %mask_sparse
            % composition_image_sparse = zeros(size(composition_image,1),size(composition_image,2),size(composition_image,3));
            % for i = 1:2
            %     composition_image_sparse(:,:,i) = sparse_main(composition_image(:,:,i),1000,0,5,100,5);
            % end
            % composition_image_sparse(:,:,3) = composition_image(:,:,3);
            composition_image_resize = zeros(finter*size(composition_image,1),finter*size(composition_image,2),size(composition_image,3));
            composition_image_resize = abs(fourierInterpolation(composition_image ,[finter,finter,1],'lateral'));
            composition_imagel = composition_image_resize;
            middle1 = composition_imagel(:,:,1);
            middle2 = composition_imagel(:,:,2);
            middle2(middle1 > middle2) = middle2(middle1 > middle2)/par1;
            middle1(middle2 > middle1) = middle1(middle2 > middle1)/par1;
            composition_image_resize (:,:,1) =  middle1;
            composition_image_resize (:,:,2) =  middle2;
            order = 2;
            SACD_2channels;
            % end
        end
    end
end

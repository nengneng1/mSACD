for nn = 1:30
    ch1_(:,:,nn) = ch1(nn,:);
    ch2_(:,:,nn) = ch2(nn,:);
    ch3_(:,:,nn) = ch3(nn,:);
end
for ga = [0.3,0.5,0.7,1]
    for mul =[0];
        lag = 1;
        % for o = [0.3,0.4,0.5,0.6,0.7,1]
        FWHM2 = 3;
        iter2 = 5;%iteration times of post decovnolution
        Ipsf2 = generate_psfv0(FWHM2*finter);

        Ipsf_mito = generate_psfv0(3*finter);
        sparse_on = 1;
        gamma2 = ga;

        normax = 100;
        %% 2 sparse param
        fidelity0 = 200;%保真度，越低则越平滑
        fidelity_z0 = 0; %z轴保真度，各项同性默认为1
        sparsity0 = 0.1; %稀疏度
        back_sparse =0;%背景参数：当设置为0时，不滤除背景；当设置为大于1的数时，越接近1滤除效果越明显


        %% SACD cumulant
        % %
        if TMIon == 1
            maskimg1 = (composition_image_resize(:,:,1).*ch1_.*stack4SACD_resize.^(gamma2));
            maskimg2 = (composition_image_resize(:,:,2).*ch2_.*stack4SACD_resize.^(gamma2));
            maskimg3 = (composition_image_resize(:,:,3).*ch2_.*stack4SACD_resize.^(gamma2));
            maskimg1 = maskimg1./max(max(maskimg1(:)));
            maskimg2 = maskimg2./max(max(maskimg2(:)));
            maskimg3 = maskimg3./max(max(maskimg3(:)));

            imwritestack(double(maskimg1.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2),num2str(FWHM4gauss3)],'_order',num2str(gamma_1),'_',num2str(gamma2),'_mask1.tif']);
            imwritestack(double(maskimg2.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2),num2str(FWHM4gauss3)],'_order',num2str(gamma_1),'_',num2str(gamma2),'_mask2.tif']);
            imwritestack(double(maskimg3.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2),num2str(FWHM4gauss3)],'_order',num2str(gamma_1),'_',num2str(gamma2),'_mask3.tif']);


            stacksub1 = abs(maskimg1 - mul * mean(maskimg1,3));%获得涨落信号
            stacksub2 = abs(maskimg2 - mul * mean(maskimg2,3));%获得涨落信号
            stacksub3 = abs(maskimg3 - mul * mean(maskimg3,3));%获得涨落信号
            stacksub1 =  stacksub1(:,:,1:10);
            stacksub2 =  stacksub2(:,:,1:20);
            stacksub3 =  stacksub3(:,:,11:30);

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
            %% 每个通道单独调节解卷积参数
            %1
            sparse1 = sparse_main(cumimg1.^0.5,80,0,1,100,back_sparse);
            SACD1 = double(abs(deconvlucy(sparse1.^2, Ipsf2.^order, 3))); %the second RL deconvolution
            SACD1 = double(SACD1./max(SACD1(:)));
            SACD1 = SACD1.^0.5;
            %2
            sparse2 = sparse_main(cumimg2.^0.5,200,0,2,100,back_sparse);
            SACD2 = double(abs(deconvlucy(sparse2.^2, Ipsf_mito.^order, 3))); %the second RL deconvolution
            SACD2 = double(SACD2./max(SACD2(:)));
            SACD2 = SACD2.^0.5;
            %3
            sparse3 = sparse_main(cumimg3.^0.25,100,0,0.1,100,3);
            SACD3 = double(abs(deconvlucy(sparse3.^4, Ipsf2.^order, 8))); %the second RL deconvolution
            SACD3 = double(SACD3./max(SACD3(:)));
            SACD3 = SACD3.^0.5;
            %%
            [x,y] = size(SACD1)
            SACD2 = imresize(SACD2,[x,y]);
            SACD3 = imresize(SACD3,[x,y]);
            SACD_composi(:,:,1) =  SACD1;
            SACD_composi(:,:,2) =  SACD2;
            SACD_composi(:,:,3) =  SACD3;
            imwritestack(double(SACD_composi.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2),num2str(FWHM4gauss3)],'_order',num2str(gamma_1),'_',num2str(gamma2),'_TMISACD.tif']);
        end
        % clear SACD_composi
        MIP = mean(imgstack(:,:,:),3);
        % WF = imresize(WF,[size(stack4SACD,1),size(stack4SACD,2)]);

        MIP = double(MIP./max(MIP(:)));
        imwritestack(double(abs(MIP.*65535)),[imgfolder,'/','back',[num2str(ifsub),'_',num2str(backg)],'_base',num2str(baseline),'_WF.tif']);
    end
end
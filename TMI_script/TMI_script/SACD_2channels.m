for nn = 1:20
    ch1_(:,:,nn) = ch1(nn,:);
    ch2_(:,:,nn) = ch2(nn,:);
end
for ga = [0.2,0.3,0.5,0.7,1,1.5]
    for mul =[0];
        lag = 1;
        % for o = [0.3,0.4,0.5,0.6,0.7,1]
        FWHM2 = 3;
        iter2 = 5;%iteration times of post decovnolution
        Ipsf2 = generate_psfv0(FWHM2*finter);
        gamma2 = ga;

        %% SACD cumulant
        % %
        if TMIon == 1
            maskimg1 = (composition_image_resize(:,:,1).*ch1_.*stack4SACD_resize.^(gamma2));
            maskimg2 = (composition_image_resize(:,:,2).*ch2_.*stack4SACD_resize.^(gamma2));
            maskimg1 = maskimg1./max(max(maskimg1(:)));
            maskimg2 = maskimg2./max(max(maskimg2(:)));
            imwritestack(double(maskimg1.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2)],'_order',num2str(gamma_1),'_',num2str(gamma2),'_mask1.tif']);
            imwritestack(double(maskimg2.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2)],'_order',num2str(gamma_1),'_',num2str(gamma2),'_mask2.tif']);
            stacksub1 = abs(maskimg1 - mul * mean(maskimg1,3));%获得涨落信号
            stacksub2 = abs(maskimg2 - mul * mean(maskimg2,3));%获得涨落信号
            stacksub1 =  stacksub1(:,:,1:20);
            stacksub2 =  stacksub2(:,:,1:20);

            Nx=size(stacksub1,1);
            Ny=size(stacksub1,2);
            cumimg1=zeros(Nx,Ny);
            cumimg2=zeros(Nx,Ny);

            cumimg1(2:Nx-1,2:Ny-1) = (mean(stacksub1(1:Nx-2,2:Ny-1,:).*stacksub1(3:Nx,2:Ny-1,:),3) ...
                + mean(stacksub1(2:Nx-1,1:Ny-2,:).*stacksub1(2:Nx-1,3:Ny,:),3))./2;
            cumimg2(2:Nx-1,2:Ny-1) = (mean(stacksub2(1:Nx-2,2:Ny-1,:).*stacksub2(3:Nx,2:Ny-1,:),3) ...
                + mean(stacksub2(2:Nx-1,1:Ny-2,:).*stacksub2(2:Nx-1,3:Ny,:),3))./2;
            %% 每个通道单独调节解卷积参数
            %1
            sparse1 = sparse_main(cumimg1.^0.5,200,0,1,100,0);
            SACD1 = double(abs(deconvlucy(sparse1.^2, Ipsf2.^order, 7))); %the second RL deconvolution
            SACD1 = double(SACD1./max(SACD1(:)));
            SACD1 = SACD1.^0.5;
            %2
            sparse2 = sparse_main(cumimg2.^0.5,200,0,2,100,0);
            SACD2 = double(abs(deconvlucy(sparse2.^2, Ipsf2.^order, 7))); %the second RL deconvolution
            SACD2 = double(SACD2./max(SACD2(:)));
            SACD2 = SACD2.^0.5;
            %%
            [x,y] = size(SACD1);
            SACD2 = imresize(SACD2,[x,y]);
            SACD1 = SACD1./max(SACD1(:));
            SACD2 = SACD2./max(SACD2(:));
            ga_ = ga*10;
            eval(['SACD_ch1_',num2str(ga_),'(:,:,f) =  SACD1;']);
            eval(['SACD_ch2_',num2str(ga_),'(:,:,f) =  SACD2;']);
        end

    end
end
imwritestack(double(SACD_ch1_2.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2)],'_order',num2str(gamma_1),'_',num2str(0.2),'_TMISACD_ch1.tif'])
imwritestack(double(SACD_ch2_2.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2)],'_order',num2str(gamma_1),'_',num2str(0.2),'_TMISACD_ch2.tif'])
imwritestack(double(SACD_ch1_3.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2)],'_order',num2str(gamma_1),'_',num2str(0.3),'_TMISACD_ch1.tif'])
imwritestack(double(SACD_ch2_3.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2)],'_order',num2str(gamma_1),'_',num2str(0.3),'_TMISACD_ch2.tif'])
imwritestack(double(SACD_ch1_5.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2)],'_order',num2str(gamma_1),'_',num2str(0.5),'_TMISACD_ch1.tif'])
imwritestack(double(SACD_ch2_5.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2)],'_order',num2str(gamma_1),'_',num2str(0.5),'_TMISACD_ch2.tif'])
imwritestack(double(SACD_ch1_7.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2)],'_order',num2str(gamma_1),'_',num2str(0.7),'_TMISACD_ch1.tif'])
imwritestack(double(SACD_ch2_7.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2)],'_order',num2str(gamma_1),'_',num2str(0.7),'_TMISACD_ch2.tif'])
imwritestack(double(SACD_ch1_10.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2)],'_order',num2str(gamma_1),'_',num2str(1),'_TMISACD_ch1.tif'])
imwritestack(double(SACD_ch2_10.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2)],'_order',num2str(gamma_1),'_',num2str(1),'_TMISACD_ch2.tif'])
imwritestack(double(SACD_ch1_15.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2)],'_order',num2str(gamma_1),'_',num2str(1.5),'_TMISACD_ch1.tif'])
imwritestack(double(SACD_ch2_15.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2)],'_order',num2str(gamma_1),'_',num2str(1.5),'_TMISACD_ch2.tif'])
imwritestack(double(abs(MIP.*65535)),[imgfolder,'/','back',[num2str(ifsub),'_',num2str(backg)],'_base',num2str(baseline),'_MIP.tif']);

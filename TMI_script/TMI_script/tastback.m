  for i =  1:1
                % BACK1(:,:,i) = background_estimation(maskimg1(:,:,i)./backg,1,6,'db6',3);
                % BACK1(BACK1 < 0) = 0;
                BACK2_(:,:,i) = background_estimation(maskimg2(:,:,i)./4,1,8,'db8',3);
                BACK2_(BACK2_ < 0) = 0;
                BACK3_(:,:,i) = background_estimation(maskimg3(:,:,i)./4,1,8,'db8',3);
                BACK3_(BACK3_ < 0) = 0;
  end
     imwritestack(double(abs(BACK3.*65535)),[imgfolder,'/',num,'_back3.tif']);
            imwritestack(double(abs(BACK2.*65535)),[imgfolder,'/',num,'_back2.tif']);
             maskimg2_ = maskimg2(:,:,1) - BACK2_*2;
            maskimg3_ = maskimg3(:,:,1)  - BACK3_*2;
            imwritestack(double(maskimg2_.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'_par',num2str(par1*10),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2),num2str(FWHM4gauss3)],'_mul',num2str(mul),'_lag',num2str(lag),'_dif',num2str(diff_param),'_order',num2str(gamma_1),'_',num2str(gamma2),'_frame',num2str(frame),'_mask2.tif']);
            imwritestack(double(maskimg3_.*65535),[imgfolder,'/','RLon',num2str(RL1_on),'back',[num2str(ifsub),'_',num2str(backg)],'_sigma',[num2str(sigma(1)),num2str(sigma(2)),num2str(sigma(3)),num2str(sigma2)],'_fwhm',[num2str(FWHM),num2str(FWHM2)],'_iter',[num2str(iter1),'_',num2str(iter2)],'_finter',num2str(finter),'_base',num2str(baseline),'_par',num2str(par1*10),'gaus',[num2str(FWHM4gauss1),num2str(FWHM4gauss2),num2str(FWHM4gauss3)],'_mul',num2str(mul),'_lag',num2str(lag),'_dif',num2str(diff_param),'_order',num2str(gamma_1),'_',num2str(gamma2),'_frame',num2str(frame),'_mask3.tif']);
            
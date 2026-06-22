clc;clear;close
for i = 7:4:79
    filefolder = ['result/202504122_3colors/',num2str(i),'_frame30']
    filename = 'RLon1back0_4_sigma1130_fwhm3.53_iter5_5_finter2_base0gaus222_order1_0.3_TMISACD.tif';
    file = [filefolder,'/',filename];
    stack=imreadstack(file);
    stack1 = stack(:,:,1);
    stack1 = stack1./max(stack1(:));
    stack1 = stack1.^2;
    stack2 = stack(:,:,2);
    stack2 = stack2./max(stack2(:));
    stack2 = stack2.^2;
    imwritestack(stack1.*65535,[filefolder,'/',num2str(i),'_nucleus','.tif']);
    imwritestack(stack2.*65535,[filefolder,'/',num2str(i),'_mito','.tif']);
    filename2 = 'RLon1back0_4_sigma1130_fwhm3.53_iter5_5_finter2_base0gaus222_order1_1_TMISACD.tif';
    file2 = [filefolder,'/',filename2];
    stack_=imreadstack(file2);
    stack3 = stack_(:,:,3);
    stack3 = stack3./max(stack3(:));
    stack3 = (stack3.^0.5);
    imwritestack(stack3.*65535,[filefolder,'/',num2str(i),'_MT','.tif']);
end

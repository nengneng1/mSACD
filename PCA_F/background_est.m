% The goal of this function is to estimate the background using either the
% end of decay or using the areas outside the region of interest
function [bg] = background_est(img,mode,system,n_ROI)
%     bg = [];

if(strcmp(system,'ISI'))
    
if (mode==0)
    
%     bg=squeeze(mean(img(:,:,100:105),3));
    bg=squeeze(mean(img(:,:,207:211),3));

end 

if (mode==1)
    
    %% cell masking
sample=img(1:32,1:32,:);
rep=repmat(sample,8);
bg=mean(rep,3);

end
if (mode==2)
    
 av_sig=(squeeze(mean(img(:,:,:),[1 2])));
 av_sig_ROI=av_sig(1:208,:);   
 figure()
 plot(av_sig_ROI)
[xData, yData] = prepareCurveData( [], av_sig_ROI );

% Set up fittype and options.
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.StartPoint = [0.0933148869465796 -0.042605233140757 0.0521349224234699 -0.0102149526878041];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'av_sig_ROI', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
ylabel av_sig_ROI
grid on

coeffs=coeffvalues(fitresult);
x=1:208;
func=coeffs(1).*exp(coeffs(2).*x)+coeffs(3).*exp(coeffs(4).*x);

DC=mean(func(190:208));
figure()
plot(x,av_sig_ROI-DC)
bg=DC.*ones(256);
end
end
if(strcmp(system,'EIT'))
 av_sig=(squeeze(mean(img(:,:,:),[1 2])));
 
 
  figure()
  plot(av_sig)
  
 av_sig_ROI=av_sig(1:n_ROI,:);   
 
[xData, yData] = prepareCurveData( [], av_sig_ROI );

% Set up fittype and options.
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.StartPoint = [0.0933148869465796 -0.042605233140757 0.0521349224234699 -0.0102149526878041];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'av_sig_ROI', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
ylabel av_sig_ROI
grid on

coeffs=coeffvalues(fitresult);
x=1:n_ROI;
func=coeffs(1).*exp(coeffs(2).*x)+coeffs(3).*exp(coeffs(4).*x);

DC=mean(func(n_ROI-5:n_ROI));
figure()
plot(x,av_sig_ROI-DC)
bg=DC.*ones(256);
end
end
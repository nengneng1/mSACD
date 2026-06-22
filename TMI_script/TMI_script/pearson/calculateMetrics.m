function [pearsonVal, ssimVal, psnrVal] = calculateMetrics(imSplit, imTruth)
vecSplit = imSplit(:);
vecTruth = imTruth(:);
pearsonVal = corr(vecSplit, vecTruth);
ssimVal = ssim(imSplit, imTruth,'Exponents', [1, 1, 1]);
psnrVal = psnr(imSplit, imTruth);
end
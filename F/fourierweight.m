function [output] = fourierweight(input,alpha)

baseSofiXFourier = input;
nx = size(baseSofiXFourier,1);
ny = size(baseSofiXFourier,2);
mdf = generate_psf(220/30,(nx-1)/2);

mxInf = floor((nx-1)/2);
myInf = floor((ny-1)/2);
mxSup = ceil((nx-1)/2);
mySup = ceil((ny-1)/2);

otm1 = abs(fftshift(fft2(mdf)));
otm1 = otm1/max(otm1(:));
otm2 = mConv2(otm1,otm1);
otm2 = otm2/max(otm2(:));
otm3 = interp2((-2*myInf:2:2*mySup)',-2*mxInf:2:2*mxSup,otm1,(-myInf:mySup)',-mxInf:mxSup);
otm3 = otm3/max(otm3(:));

weight = otm3./(alpha+otm2); 

for j = 1:size(input,3)
    imageSofiFourier(:,:,j) = abs(ifft2(ifftshift(fftshift(fft2(baseSofiXFourier(:,:,j))).*weight)));
end
imageSofiFourier = double(imageSofiFourier);
output = imageSofiFourier./max(max(max(imageSofiFourier)));
    









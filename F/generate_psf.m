function y=generate_psf(gama,kernelRadius)
gama=gama/(2.335);
if nargin == 1 
    kernelRadiusl = ceil(gama* sqrt(-2* log(0.0002)))+1;
    kernelRadius = min(kernelRadiusl, 300/2 - 1);
end
ii=-kernelRadius:kernelRadius;
rsf_x=1/2*(erf((ii+0.5)./(sqrt(2).*gama))-erf((ii-0.5)./(sqrt(2).*gama)));
kernel= rsf_x'* rsf_x;
y=kernel./sum(kernel(:));
end
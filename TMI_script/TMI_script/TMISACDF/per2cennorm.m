function output = per2cennorm(data, miper, maper)

if nargin < 3, maper = 100; end
if nargin < 2, miper = 0; end
data = single(data);

for i = 1 : size(data,3)
    mati = data(:,:,i);
    datamin = prctile(mati(:), miper);
    datamax = prctile(mati(:), maper);
    output(:,:,i) = (data(:,:,i) - datamin) / (datamax - datamin);
end
output(output > 1) = 1;
output(output < 0) = 0;

end
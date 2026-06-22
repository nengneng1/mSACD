function [output] = Twavelet(input)

im = input;
for k = 1:size(im,1)
    for j = 1:size(im,2)
        ob = squeeze(im(k,j,:));
[c,s] = wavedec(ob,7,'db4');
c(1:s(1)*1) = 0;
ob2 =  waverec(c,s,'db4');
im(k,j,1:size(im,3)) = ob2;
    end
end

output = im./max(max(max(im)));

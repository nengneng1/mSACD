function [output] = Xwavelet(input)
% im = fourierInterpolation(input,[2,2,1],'lateral');
im = input;
for i = 1:size(im,3)
[c,s] = wavedec2(im(:,:,i),3,'db4');
% c(end-s(8)*s(8)*3:end) = 0;
c(end-s(end-1)*s(end-1)*3:end) = 0;
ob2(:,:,i) =  waverec2(c,s,'db4');
ob = ob2;
% ob = imresize(ob2,[size(input,1),size(input,2)]);
end
% size(c)
% s(1)*s(1)*1+s(2)*s(2)*3+s(3)*s(3)*3+s(4)*s(4)*3+s(5)*s(5)*3+s(6)*s(6)*3+s(7)*s(7)*3+s(8)*s(8)*3
output = ob./max(max(max(ob)));

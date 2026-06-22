function image=gencavelaes(imglarge,x,y,r1,r2) 
%   r1，r2 分别为外半径和内半径，r1应大于r2
%  x,y>r1
c=genRing(2*r1+1,r1,r2);
imglarge(x-r1:x+r1,y-r1:y+r1)=imglarge(x-r1:x+r1,y-r1:y+r1)+c;
image=imglarge;
end
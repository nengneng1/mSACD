function [ c ] = genCircle(w,r)
%GENCIRCLE Summary of this function goes here
%   Detailed explanation goes here
%   w 是模板的大小
%   r 圆形模板的半径
[rr cc] = meshgrid(1:w);
c = sqrt((rr-floor(w/2)).^2 + (cc-floor(w/2)).^2) <= r;
end
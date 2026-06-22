function [ c ] = genRing( w,r1,r2 )
%GENRING Summary of this function goes here
%   Detailed explanation goes here
%   r1，r2 分别为外半径和内半径，r1应大于r2
c1 = genCircle(w,r1);
c2 = genCircle(w,r2);
c = c1 & ~c2;
end
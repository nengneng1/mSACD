function result = covariance(input)
input = abs(input-1*mean(input,3));
skip = size(input,3);
M11 = mean(input(:,:,1:end-1) .* input(:,:,2:end),3);
M02 = mean(input(:,:,2:end).^2,3);
M20 = mean(input(:,:,1:end-1).^2,3);
M22 = mean((input(:,:,1:end-1).^2) .* (input(:,:,1:end-1).^2),3);
result = -(skip - 2) * M11.^2 / ((skip * (skip-1))) + M02.*M20 / (skip*(skip-1)) + M22 / skip;

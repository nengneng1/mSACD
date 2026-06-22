
function [lower_bound,upper_bound] = bin_finder(bin_def,number)
lower_bound=[];
upper_bound=[];
if (numel(number)==1)
% Find which bin the random number falls into
[~, bin_index] = min(abs(bin_def - number));

% Determine the upper and lower boundaries of the bin
if bin_index == 1
    lower_bound = -1; % Lower boundary is -1
    upper_bound = bin_def(bin_index);
elseif bin_index == numel(bin_def)
    lower_bound = bin_def(bin_index - 1);
    upper_bound = 1; % Upper boundary is 1
else
    
    if (bin_def(bin_index)>number)        
    lower_bound = bin_def(bin_index - 1);
    upper_bound = bin_def(bin_index);
    else
    lower_bound = bin_def(bin_index);
    upper_bound = bin_def(bin_index+1);   
    end
end

else

for kk=1:numel(number)
    
  [~, bin_index] = min(abs(bin_def - number(kk)));

% Determine the upper and lower boundaries of the bin
if bin_index == 1
    lower_bound =[lower_bound,-1]; % Lower boundary is -1
    upper_bound =[upper_bound, bin_def(bin_index)];
elseif bin_index == numel(bin_def)
    lower_bound =[lower_bound, bin_def(bin_index - 1)];
    upper_bound = [upper_bound,1]; % Upper boundary is 1
else
    
    if (bin_def(bin_index)>number)        
    lower_bound =[lower_bound, bin_def(bin_index - 1)];
    upper_bound = [upper_bound,bin_def(bin_index)];
    else
    lower_bound = [lower_bound,bin_def(bin_index)];
    upper_bound =[upper_bound, bin_def(bin_index+1)];   
    end
end  
    
    
    
    
    
    
    
end   
    
end    
end
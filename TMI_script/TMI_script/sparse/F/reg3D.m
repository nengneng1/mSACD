function [output2] = reg3D(input1,input2)
Target = (input1);
Source = (input2);
[D,moving_reg] = imregdemons(Source, Target,'PyramidLevels',2,'AccumulatedFieldSmoothing',1.5);
% output1 = single(gather(Target));
output2 = single((moving_reg));

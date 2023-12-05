function [Fdp,r] = GetCDTf(F1,F2)
% [Fdp,r] = GetCDTf(F1,F2);
% Get Cubic DT from F1 and F2 and their ratio. 
% Just for ease in calculation

Fdp = 2*F1-F2;
r = F2/F1;

end
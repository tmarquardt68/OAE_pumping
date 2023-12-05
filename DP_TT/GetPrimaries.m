function [Set1, Set2, Set3, T] = GetPrimaries(Fdp,r);
% [Set1, Set2, Set3, T] = GetPrimaries(Fdp,r); 
% Get primary set from desired DPOAE (cubic DT) frequency and F2/F1 ratio
% Fdp : Cubic difference tone frequency. Enter desired one from which 2 others (surrounding it) 
% will be calculated in +- 1 ERB steps. 
% r: frequency ratio (e.g. 1.2, 1.22)
% SetX give: [F1 F2 Fdp] values 
% Table is also generated. Type "T" in workspace
% Fdp should preferably be multiple of 5 
% Output F1 and F2 (and CDT frequencies) are multiples of 5. 
% [exact F2/F1 ratio is usually sacrificed (but only a bit) to get this.] 

ERB = 24.7*(4.37*Fdp/1000+1); % ERB formula (Glasberg and Moore). To use approx. 1 ERB steps
Ratio = ERB/Fdp;
FSet = [round(Fdp-Ratio*Fdp) round(Fdp) round(Fdp+Ratio*Fdp)];
%      [  - 1 ERB               0            + 1 ERB        ]

for i = 1:3
F1(i) = FSet(i)/(2-r); 
F2(i) = r*F1(i);
F1(i) = F1(i)-rem(F1(i),5);
F2(i) = F2(i)-rem(F2(i),5);
end

Set1 = [F1(1) F2(1) 2*F1(1)-F2(1)]; % One way of output, into vectors ([F1 F2 Fdp])
Set2 = [F1(2) F2(2) 2*F1(2)-F2(2)]; 
Set3 = [F1(3) F2(3) 2*F1(3)-F2(3)];

CubicDT = [Set1(3);Set2(3);Set3(3)];
F1 = F1';
F2 = F2';
Rows = {'-1 ERB';'Centred' ;'+1 ERB'};
T = table(F1,F2,CubicDT,'RowNames',Rows); % Other way of output, into table

end


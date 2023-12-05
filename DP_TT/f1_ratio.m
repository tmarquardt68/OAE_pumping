function [f1_dp,f1_sf,l1_sf] = f1_ratio(f2)

r = 1.02305; 
for n=1:17, 
    f1_dp(n)=round(f2/(r.^n)/5)*5; 
    if rem(f1_dp(n),10)==0
        f1_dp(n) = f1_dp(n) + 5;
    end
    dp=2*f1_dp(n)-f2; 
    dp2= 3*f1_dp(n)-2*f2; 
    
    f1_sf(n) = f2 - dp;
    if rem(f1_sf(n),10)~=0
        f1_sf(n) = f1_sf(n) + 5;
    end    
    [n r.^n f1_dp(n) dp/1000 dp2/1000];
end, 
f1_dp;
f1_sf = f1_sf(1:13);
l1_sf = 90 - log2(f1_sf)*6 +log2(120)*6;
l1_sf = round(l1_sf(1:13));
% semilogx(f1_sf,l1_sf,'x-'), grid

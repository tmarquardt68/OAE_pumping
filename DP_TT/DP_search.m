function DP_search(f2_f1)

f1 = f2_f1*5;
if rem(f1,10)==0
    f1 = f1+5;
end
f1 = [f1-200 f1 f1+200];
f2 = round(f1*1.2/5)*5;
f2_f1 = f2-f1

f2 = repmat(f2,1,4)
f1 = repmat(f1,1,4)



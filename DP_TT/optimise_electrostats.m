function [harm_1,harm_2,harm_3,harm_4,harm_5, H] = optimise_electrostats

fs = 48000; 
tone = sin(linspace(0,1000*2*pi*(1-1/fs), fs));
load('C:\Users\Admin\Documents\MATLAB\OAE_gerbil\DP_TT\DP_TT_calibration_data\00000__NONE__.mat','H_mic')

V_pol = 198:202; % around 200V
V_offset = 212.1:0.2:213.5; % around 212.7V
H = NaN*ones(length(V_pol),length(V_offset),fs);
harm_1 = NaN*ones(length(V_pol),length(V_offset));
harm_2 = harm_1;
harm_3 = harm_1;
harm_4 = harm_1;
harm_5 = harm_1;

for q = 1:length(V_pol)
    for q2 = 1:length(V_offset)
        tone_out = ((V_pol(q)*tone+V_offset(q2)).^0.5-sqrt(V_offset(q2)))/sqrt(V_pol(q));
        data_in = data_acquisition(tone_out', 1, 1, 1);
        %H(q,q2,:) = 20*log10(abs(fft(data_in)./H_mic));
        H(q,q2,:) = 20*log10(abs(fft(data_in)));
        harm_1(q,q2) = H(1001);
        harm_2(q,q2) = H(2001);
        harm_3(q,q2) = H(3001);
        harm_4(q,q2) = H(4001);
        harm_5(q,q2) = H(5001);
        stem(2:10000,squeeze(H(q,q2,2:10000)+20)),ylim([0 100]),grid,drawnow
        [q q2]
    end
end

% V_pol = 190:198; % around 200V
% V_offset = 212.5:0.1:212.9; % around 212.7V
% 
%   -11.1865  -12.2872  -13.2582  -11.8962  -10.8185
%   -12.3554  -12.9381  -12.8854  -11.6276  -11.0797
%   -12.8580  -11.3345  -11.8346  -10.9743  -10.1535
%   -13.1206  -12.3269  -12.2123  -10.6046  -10.8677
%   -11.8149  -12.3678  -10.7555  -12.6198  -11.6063
% 
%     -21.9231  -21.2824  -26.1571  -23.9431  -23.2124
%   -21.5947  -20.6204  -19.3127  -41.9766  -27.1235
%   -22.9352  -21.9354  -23.9952  -23.1039  -24.4646
%   -18.5839  -23.1408  -20.3425  -22.1749  -21.2327
%   -21.3397  -20.5210  -20.2714  -23.7477  -31.2717


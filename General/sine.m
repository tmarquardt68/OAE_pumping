% SINE returns a matrix with rows of sine waves Amplitude is set to +/- 1.
% 
%   wave = SINE(sample_rate,frequencies):
%
%   $sample_rate:  waveforms will calculated for a given sampling frequency.
%   $frequencies:  vector containing frequency values:
%                	- 0 Hz outputs a row of zeros
%                	- negativ values give 180 degree phase inversion
%
% The length of the matrix depends on the greatest common divisor of
% the frequencies $f1...$f8. So the rows contain an integral number of 
% sine cicles (usable in a loop).
%
% SINE is part of the PC_AUDIO toolbox (by Torsten Marquardt)

function wave = sine(sample_rate,frequencies)

error(nargchk(2,2,nargin)),
n_freg = length(frequencies);

l = sample_rate;
for (n = 1:n_freg)
    l = gcd(l,abs(frequencies(n)));
end
l = sample_rate/l;
for (n = 1:n_freg)
    wave(n,:) = linspace(0,frequencies(n)*...
        l/sample_rate*2*pi*(1-1/l), l);
end
wave = sin(wave);

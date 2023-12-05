% [buffer_in, n_used] = DATA_ACQUISITION ( ...
%							wave_out, out_chs, in_chs, n_avg, use_crit)
%
% The $wave data are played $n_avg+2 times and first and the last loop 
% is discarded and therefore not averaged.
%
% Part of the OAE toolbox 
% Copyright (C) 2008 Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).


function [buffer_in, n_used] = data_acquisition( ...
	wave_out, out_chs, in_chs, n_avg, use_crit)

global SAMPLE_RATE,

if isempty(SAMPLE_RATE),
	SAMPLE_RATE = 48000;
end,

le = length(wave_out);
h_sound_device = sound_io('open_device',SAMPLE_RATE, out_chs, in_chs);
wave_in = sound_io('playrecord',h_sound_device,...
    repmat(wave_out,(n_avg+3),1),(n_avg+3)*le/SAMPLE_RATE);
sound_io('close_device',h_sound_device);
wave_in = reshape(wave_in(2*le+1:(n_avg+2)*le,:), le, n_avg,[]);
buffer_in = double(squeeze(sum(wave_in, 2)./n_avg));
n_used = 0;



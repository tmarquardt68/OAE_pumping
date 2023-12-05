function varargout = sound_io(varargin);
%
% Sound input output routine using ASIO of PsychPort_Audio.mex
% (in PsychToolbox)
% 
% v1.0 3-08-12 (UCL Ear Institute)
%
%
% Part of the Stimulus Toolbox 
% Copyright (C) Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).
%
% TO DO
% - playrecord(): return last the appropriate amountToAllocateSecs bytes of
%   the 1.5*amountToAllocateSecs (compensatinf for dealy)

% FEVAL switchyard
if (nargout),
    [varargout{1:nargout}] = feval(varargin{1},varargin{2:end});
else,
    feval(varargin{1},varargin{2:end}); 
end,
%==========================================================================
%==========================================================================
%
%==========================================================================
function h_sound_device = open_device(sampleRate, nOutputCh, nInputCh)

% Specify here a certain sound device by hand. 
% If left empty, the first ASIO device found will be used.
deviceid = [];

buffersize = 1024;     % Pointless to set this. Auto-selected to be optimal.
reqlatencyclass = 1;
% selectchannels = [0:63;[0:63]+4]; % for channel mapping.!!! Currently rec shifted!! 
selectchannels = [0:63;0:63]; % for channel mapping.!!! Currently rec shifted!! 

h_sound_device = -1;

if nOutputCh > 0 && nInputCh > 0 % playrec
    selectchannels = selectchannels(:,1:max(nInputCh,nOutputCh));
    channels = [nOutputCh nInputCh];
    mode = 3;
elseif nOutputCh > 0 % playback only
    selectchannels = selectchannels(1,1:nOutputCh);
    channels = nOutputCh;
    mode = 1;
elseif nInputCh > 0 % record only
    selectchannels = selectchannels(2,1:nInputCh);
    channels = nInputCh;
    mode = 2
else
    error(['nOutputCh =  ' num2str(nOutputCh) ...
        ', nInputCh =  ' num2str(nInputCh)  '!'])
end

% Link PsychPortAudio against the portaudio DLL
olddir = pwd;
drv_path = fileparts(which('portaudio_x86.dll'));
cd(drv_path);
d = PsychPortAudio('GetDevices', 3);
cd(olddir);

if ~isempty(d)
    % And some more commercials as required by the license...
    %fprintf('Using "ASIO Interface Technology by Steinberg Media Technologies GmbH"\n\n');
    if isempty(deviceid)
        deviceid = d(1);
    end
    % selectchannels = []; % channel mapping not required
    h_sound_device = PsychPortAudio('Open',deviceid.DeviceIndex,mode,...
        reqlatencyclass,sampleRate,channels,buffersize,[],selectchannels,16);
    
else
    errordlg('No ASIO sound devices found!')
end

%%==========================================================================
function play(h_sound_device,bufferdata)

PsychPortAudio('FillBuffer',h_sound_device,bufferdata');
PsychPortAudio('Start', h_sound_device);

%==========================================================================
function err = record(fs)

if playrec('getSampleRate')~= fs,
    error('Sample_rate mismatch!')
end
error('Not supported yet!')

%==========================================================================
function recBuffer = playrecord(h_sound_device,bufferdata,amountToAllocateSecs)

PsychPortAudio('FillBuffer',h_sound_device,bufferdata');
PsychPortAudio('GetAudioData', h_sound_device, 1.5*amountToAllocateSecs);
PsychPortAudio('Start', h_sound_device);
recBuffer = PsychPortAudio('GetAudioData',h_sound_device, [],amountToAllocateSecs)';
%==========================================================================
function err = clear_up(h_sound_device)
% Kept for compatibility with version soiund_io_playrec.m (where this
% function is necessary to flush buffers(?))

err = 0;

%==========================================================================
function err = close_device(h_sound_device)

err = 0;
if ~exist('h_sound_device','var')|| isempty(h_sound_device)
    PsychPortAudio('Close') %close all open devices
else 
    PsychPortAudio('Stop', h_sound_device);
    PsychPortAudio('Close', h_sound_device)
end

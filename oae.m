% OAE.M :  starts the OAE-System.
% Please change this values according to your specific system

% Part of the OAE toolbox 
% Copyright (C) 2008 Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).


global SAMPLE_RATE BITS_PER_SAMPLE OAE_PATH MIXER

% Directory where the OAE-Toolbox is situated:
OAE_PATH = 'C:\Users\Admin\Documents\MATLAB\OAE_gerbil';

% some soundcard specific informations:
SAMPLE_RATE = 48000;
BITS_PER_SAMPLE = 24;

% What is the mixer application you want to use to see what happens at 
% the sound port? (Usefull only if it shows level meter)
% empty string => no mixer 
MIXER = ''; %'C:\Program Files\MOTU\Audio\Apps\MOTU CueMix Console.exe';

addpath([OAE_PATH,'\General\'],'-begin'),

% Here it starts
oae_subject_browser('ini'),

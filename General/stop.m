% Part of the OAE toolbox 
% Copyright (C) 2008 Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).

function stop(action)

switch action,
    case 'ini',
        stop_fig,
    case 'stop',
        set(findobj('Tag','continue'),'Tag','stop'),
        set(findobj('Name','Stop'),'Name','Series will stop!'),
end,
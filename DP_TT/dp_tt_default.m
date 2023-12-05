% DP_TT_DEFAULT(action):
% lets you browse chose and delete default parameter.
%
% Part of the OAE toolbox 
% Copyright (C) 2008 Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).

function dp_tt_deault(action)

global OAE_PATH 

switch action
    case 'ini'
        h = findobj('Name','DP_TT Default Browser');
        if(isempty(h)),
            dp_tt_default_fig,
            h = findobj(gcf,'Tag','defaults');
            set(h,'Value',1),
            dp_tt_default('fill')
        else,
            figure(h)
        end,
        
        
    case 'fill'
        % fill boxes with defaults
        if (~exist([OAE_PATH,'\DP_TT\DP_TT_default_data']))
            eval(['!mkdir ', OAE_PATH,'\DP_TT\DP_TT_default_data'])
        end
        str = cellstr(my_dir([OAE_PATH '\DP_TT\DP_TT_default_data'],'mat'));
        h = findobj(gcf,'Tag','defaults');
        set(h,'String',str),
        val = get(h, 'Value');
        if (isempty(str)),
            set(findobj(gcf,'Tag','F1'),'String',''),
            set(findobj(gcf,'Tag','F2'),'String',''),
            set(findobj(gcf,'Tag','F_TT'),'String',''),
            set(findobj(gcf,'Tag','L1'),'String',''),
            set(findobj(gcf,'Tag','L2'),'String',''),
            set(findobj(gcf,'Tag','L_TT'),'String',''),
            set(findobj(gcf,'Tag','nested_mode'),'String',''),
            set(findobj(gcf,'Tag','n_avg'),'String',''),
            set(findobj(gcf,'Tag','created'),'String',''),
            set(findobj(gcf,'Tag','purpose'),'String',''),
            set(findobj(gcf,'Tag','rename'),'String',''),
        else
            item = str{val};
            try,
               load([OAE_PATH,'\DP_TT\DP_TT_default_data\',item]),
               set(findobj(gcf,'Tag','F1'),'String',F1),
               set(findobj(gcf,'Tag','F2'),'String',F2),
               set(findobj(gcf,'Tag','F_TT'),'String',F_TT),
               set(findobj(gcf,'Tag','L1'),'String',L1),
               set(findobj(gcf,'Tag','L2'),'String',L2),
               set(findobj(gcf,'Tag','L_TT'),'String',L_TT),
               set(findobj(gcf,'Tag','nested_mode'),'String',nested_mode),
               set(findobj(gcf,'Tag','n_avg'),'String',t_avg),
               set(findobj(gcf,'Tag','created'),'String',created),
               set(findobj(gcf,'Tag','purpose'),'String',comment),
               set(findobj(gcf,'Tag','rename'),'String',str{val}),
            catch,
            end,
         end
         
        
    case 'load'
        % fill boxes of dp_TT_measurement window with default parameter
        h = findobj(gcbf,'Tag','defaults');
        str = get(h, 'String');
        if (~isempty(str)), 
            val = get(h, 'Value');
            item = str{val};                
            dp_tt_measurement('load_defaults', item), 
        else
            dp_tt_measurement('load_defaults', 'none'),          
        end,
        close(gcf),
        
        
    case 'delete'
        % delete default parameter set choosen in listbox
        h = findobj(gcbf,'Tag','defaults');
        str = get(h, 'String');
        val = get(h, 'Value');
        item = str{val};                
        if (isempty(item)), return, end,
        if (strcmp(questdlg(strvcat('Are you sure you want delete : ',...
            item),'WARNING - DP_TT Default Browser:','No'),'Yes')),
            delete ([OAE_PATH,'\DP_TT\DP_TT_default_data\', item])
			dp_tt_default('fill'),
        end 
        
        
    case 'rename'
        h = findobj(gcbf,'Tag','defaults');
        val = get(h,'Value');
        str = get(h,'String');
        load([OAE_PATH,'\DP_TT\DP_TT_default_data\',str{val}]),
        default = get(findobj(gcbf,'Tag','rename'),'String');
        save([OAE_PATH,'\DP_TT\DP_TT_default_data\',default],...
            'F1', 'F2', 'F_TT','L1', 'L2', 'L_TT', 'nested_mode',...
            't_avg','created','comment','default','min_freq')
        delete([OAE_PATH,'\DP_TT\DP_TT_default_data\',str{val}]),
        
    case 'comment',
        h = findobj(gcbf,'Tag','defaults');
        str = get(h,'String');
        comment = get(findobj(gcbf,'Tag','purpose'),'String');
        save([OAE_PATH,'\DP_TT\DP_TT_default_data\',...
            str{get(h, 'Value')}],'comment','-append'),

    otherwise
        error('Unkown action (case)!');       
end %switch action

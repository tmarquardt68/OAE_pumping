% DP_TT_BROWSER(action,subject):
%
% The 'subject' will be saved as part of  the figure name
% The series file names are saved in the figure property 'User'
%
%
% Part of the OAE toolbox 
% Copyright (C) 2008 Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).


function dp_tt_browser(action,subject)

global OAE_PATH MIN_FREQ
switch action
    case 'ini'
        dp_tt_browser_fig,
        set(gcf,'Name',['DP_TT Browser:  ',subject]),
        dp_tt_browser('fill', subject),
        % That variable defines the spacing of the spectral lines and therefore
        % the length of the fft-window:
        MIN_FREQ = 1;
        
        
    case 'fill', % fill listbox
        h_fig = findobj('Name',['DP_TT Browser:  ',subject]);
        h = findobj(h_fig,'Tag','series');
        files = my_dir([OAE_PATH '\Subjects\', subject],'mat');
        str = ''; 
        if (isempty(files)),
            val  = 0;
        else,
            %fill cell array 'str' with date and default_name
            [n m] = size(files);
            for (n = 1:n),
                comment = ' ';
                load([OAE_PATH '\Subjects\',subject,...
                        '\',files(n,:)],'comment'),
                if(isempty(comment)),
                    str{n} = [files(n,6:20),'  '];
                else,
                    str{n} = [files(n,6:20),'  ',comment(1,:)];
                end,
            end,
            val=length(str);
        end,        
        set(h_fig,'User',files);
        set(h,'Value',val),
        set(h,'String',str),
        dp_tt_browser('fill_text',subject),
        
    
    case 'fill_text'
        if (exist('subject'))'
            h_fig = findobj('Name',['DP_TT Browser:  ',subject]);
        else,
            h_fig = gcbf;
            name = get(gcbf,'Name');
            subject = name(17:length(name));
        end,
        h = findobj(h_fig,'Tag','series');
        str = (get(h,'String'));
        val = (get(h,'Value'));
        files = get(h_fig,'User');
		if (~isempty(files)),
			load([OAE_PATH,'\Subjects\',subject,'\',files(val,:)], ...
				't_avg','created','comment','l1','ear_l'),
			set(findobj(h_fig,'Tag','created'),'String', created),
			set(findobj(h_fig,'Tag','t_avg'),'String',num2str(t_avg)),
			set(findobj(h_fig,'Tag','n_meas'),'String',num2str(length(l1(~isnan(l1))))),
			set(findobj(h_fig,'Tag','comment'),'String',comment),
			if(exist('ear_l'))
				if(ear_l),
					set(findobj(h_fig,'Tag','ear'),'String','left'),
				else,
					set(findobj(h_fig,'Tag','ear'),'String','right'),
				end,
			else,
				set(findobj(h_fig,'Tag','ear'),'String',''),
			end,
		end,

        
    case 'new'
        files = get(gcbf, 'User');
        if (isempty(files)),
            series = '';
        else
            series = files(get(findobj(gcbf,'Tag','series'), 'Value'),:);
        end
        name = get(gcbf,'Name');
        dp_tt_measurement('ini', name(17:length(name)),series),


    case 'view'
        files = get(gcbf, 'User');
        if (isempty(files)), return, end,
        series = files(get(findobj(gcbf,'Tag','series'), 'Value'),:);
        name = get(gcbf,'Name');
        h = findobj('Name',['DP_TT_Viewer: ',name(17:length(name)),' , '...
            series(6:20)]);
        if (isempty(h))
            dp_tt_viewer('ini',name(17:length(name)),series),
        else,
            figure(h)
        end,
        
        
    case 'delete',
        h = findobj(gcbf,'Tag','series');
        str = get(h, 'String');
        val = get(h, 'Value');
        if (isempty(str)), return, end,
        files = get(gcbf, 'User');
        if (strcmp(questdlg(strvcat('Are you sure you want to delete :',...
                str{val}),'WARNING - DP_TT Browser:'),'Yes')),
            name = get(gcbf,'Name');
            
            % close all open viewer windows of the series
            h_viewer = findobj('Name',['DP_TT_Viewer: ',...
                    name(17:length(name)),' , ',files(val,6:20)]);
            if (~isempty(h_viewer)),
                close(h_viewer),
            end,
            filenames = files(val,:);
            % delete file
            delete([OAE_PATH,'\Subjects\',name(17:length(name)),'\',...
                    filenames(1:20),'*.*']),
            
            % update listbox
            set(h,'Value',1);
            dp_tt_browser('fill',name(17:length(name))),
            set(h, 'Value',length(get(h,'String'))),
        end
       
         
    case 'comment'
        h = findobj(gcbf,'Tag','series');
        comment = get(findobj(gcbf,'Tag','comment'),'String');
        val = get(h,'Value');
        str = get(h,'String');
        str{val} = [str{val}(1:17),comment(1,:)];
        set(h,'String',str);
        name = get(gcbf,'Name');
        files = get(gcbf,'User');
        save([OAE_PATH,'\Subjects\',name(17:length(name)),'\',...
                files(val,:)],'comment','-append'),
        % update listbox
        dp_tt_browser('fill',name(17:length(name))),
            
    case 'probe_fit'
        h_fig = gcbf;
        val = get(findobj(gcbf,'Tag','series'),'Value');
        files = get(gcbf,'User');
        if (isempty(files)),
            warndlg('Probe fit data not available (first measurement).'...
                , 'dp_tt_browser'),
            return,
        end,
        name = get(gcbf,'Name');
        subject = name(17:length(name));
        load([OAE_PATH,'\Subjects\',subject,'\',...
                files(val,:)],'HdB_ear_channel','min_freq'),
        HdB_ear_channel = HdB_ear_channel((round(15/min_freq)+1):...
            (round(10000/min_freq)+1));

        
        % open figure
        h = findobj('Tag','magnitude_fig');
        if (isempty(h)),
            magnitude_fig,
            h = findobj('Tag','magnitude_fig');
        end,
        figure(h),
        co = get(gca,'ColorOrder');
        set(gca,'ColorOrder',[co(2:length(co),:); co(1,:)]),
        semilogx(linspace(round(15/min_freq)*min_freq,...
            round(10000/min_freq)*min_freq,length(HdB_ear_channel)-1),...
            HdB_ear_channel(1:end-1)),
        hold on, grid on, zoom on,
        
        
    case 'help ',
        
end % switch(action)

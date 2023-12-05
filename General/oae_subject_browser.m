% OAE_SUBJECT_BROWSER(action):
%
% Part of the OAE toolbox
% Copyright (C) 2008 Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).

function oae_subject_browser(action, name)

global OAE_PATH

switch action
	case 'ini'
		oae_subject_browser_fig,
		oae_subject_browser('fill'),


	case 'fill'
		% fill boxes with subject data
		if (~exist([OAE_PATH,'\Subjects']))
			eval(['!mkdir ', OAE_PATH,'\Subjects'])
		end
		files = my_dir([OAE_PATH '\Subjects']);

		if (length(files)==0),
			set(findobj(gcf,'Tag','born'),'String',''),
			set(findobj(gcf,'Tag','phone'),'String',''),
			set(findobj(gcf,'Tag','address'),'String',''),
			set(findobj(gcf,'Tag','comment'),'String',''),
			set(findobj(gcf,'Tag','n_audiogram'),'String',''),
			set(findobj(gcf,'Tag','n_dp_tt'),'String',''),
			set(findobj(gcf,'Tag','n_lf_ag'),'String',''),
			set(findobj(gcf,'Tag','dpoae'),'String',''),
			set(findobj(gcf,'Tag','teoae'),'String',''),
			return;
		end,

		str = cellstr(files);
		h_subjects = findobj(gcf,'Tag','subjects');

		val = get(h_subjects,'Value');
		set(h_subjects,'String',str(1:end)),
		str = get(h_subjects,'String');

		if exist('name'), % if subject specified in function call
			for val = 1: length(str)
				if strcmp(deblank(str{val}),deblank(name)),
					set(h_subjects,'Value',val),
				end,
			end,
		end
		
		item = str{val};
% 		if ~exist([OAE_PATH,'\Subjects\',item,'\Personal_data.mat']),
% 			born =''; phone =''; address =''; comment ='';
% 			save([OAE_PATH,'\Subjects\',item,'\Personal_data'],...
% 				'born','phone','address','comment'),
% 		end
%		load([OAE_PATH,'\Subjects\',item,'\personal_data']),
% 		set(findobj(gcf,'Tag','born'),'String', born),
% 		set(findobj(gcf,'Tag','phone'),'String',phone),
% 		set(findobj(gcf,'Tag','address'),'String',address),
%		set(findobj(gcf,'Tag','comment'),'String',comment),
		mat_files = getfield(what([OAE_PATH,'\Subjects\',item]),'mat');
		n_audiogram=0; n_dp_tt=0; n_lf_ag=0; n_dpoae=0; n_teoae=0;
		for (n=1:length(mat_files)),
			file_name = mat_files{n};
			if (file_name(1:4) == 'aud_'),
				n_audiogram = n_audiogram+1;
			elseif (file_name(1:5) == 'dp_tt'),
				n_dp_tt = n_dp_tt+1;
			elseif (file_name(1:5) == 'lf_ag'),
				n_lf_ag = n_lf_ag+1;
			elseif (file_name(1:5) == 'dpoae'),
				n_dpoae = n_dpoae+1;
			elseif (file_name(1:5) == 'teoae'),
				n_teoae = n_teoae+1;
			end,
		end,
		set(h_subjects,'String',str),
		set(findobj(gcf,'Tag','n_audiogram'),'String',n_audiogram),
		set(findobj(gcf,'Tag','n_dp_tt'),'String',n_dp_tt),
		set(findobj(gcf,'Tag','n_lf_ag'),'String',n_lf_ag),
		set(findobj(gcf,'Tag','n_dpoae'),'String',n_dpoae),
		set(findobj(gcf,'Tag','n_teoae'),'String',n_teoae),


	case 'browse'
		h = findobj(gcbf,'Tag','subjects');
		str = get(h, 'String');
		if (isempty(str)), return, end,
		subject = str{get(h, 'Value')};

		h = findobj(gcbf,'Style','radiobutton');
		for (n=1:length(h)),
			if (get(h(n),'Value') == 1),
				str =get(h(n),'String');
				eval('cd([OAE_PATH,''\'',str(1:5),])','err=1'),
				if(exist('err')),
					msgbox ('Choosen method not correctly installed!',...
						'ERROR: OAE-Subject Browser','error'),
					return,
				end,
				path([OAE_PATH,'\', str(1:5)],path),
				h_browser_fig = findobj('Name',[str(1:5),' Browser:  '...
					,subject]);
				if (isempty(h_browser_fig))
					str = [lower(str(1:5)),'_browser(''ini'',', '''' subject '''',')'];
					eval(str),
				else,
					figure(h_browser_fig)
				end,
			end,
		end

	case 'new'
		subject = inputdlg(strvcat('Please enter new patien name:',...
			'    surname_prename (seperated by underline!)'));
		if(isempty(subject)), return, end,
		subject=subject{1};
		if (isempty(subject)), return, end,
		if (exist([OAE_PATH,'\Subjects\',subject]))
			msgbox('Subject already in data bank!',...
				'ERROR OAE Subject Browser:','error');
			return
		else
			eval([subject,'=1;'],'err =1;')
			if (exist('err')),
				msgbox('Invalid file name! Use underlines instead of spaces!'...
					,'ERROR OAE Subject Browser:','error');
				return,
			end,
		end
		eval(['mkdir(''',OAE_PATH,'\Subjects\',subject ''')']),
% 		born =''; phone =''; address =''; comment ='';
% 		save([OAE_PATH,'\Subjects\',subject,'\Personal_data'],...
% 			'born','phone','address','comment'),
		oae_subject_browser('fill', subject)


	case 'delete'
		h = findobj(gcbf,'Tag','subjects');
		str = get(h, 'String');
		if (isempty(str)), return, end,
		item = str{get(h, 'Value')};
		if (strcmp(questdlg(strvcat('Sorry - Not implemented yet.',...
				item),'WARNING - OAE Subject Browser:'),'Yes')),
		end


	case 'save'
		h = findobj(gcf,'Tag','subjects');
		str = get(h,'String');
		item = str{get(h, 'Value')};
		born =    get(findobj(gcbf,'Tag','born'),'String');
		phone =   get(findobj(gcbf,'Tag','phone'),'String');
		address = get(findobj(gcbf,'Tag','address'),'String');
		comment = get(findobj(gcbf,'Tag','comment'),'String');
		save([OAE_PATH,'\Subjects\',item,'\personal_data'],...
			'born','phone','address','comment'),


	case 'radiobutton_press'
		set(findobj(gcbf,'Style','radiobutton'),'Value',0)
		set(gcbo,'Value',1)


	case 'help'

	otherwise
		error('Unkown action (case)!');
end % switch(action)


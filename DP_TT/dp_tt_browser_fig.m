function dp_tt_browser_fig()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.

load dp_tt_browser_fig                 

a = figure('Color',[0.75 0.75 0.75], ...
	'Colormap',mat0, ...
	'MenuBar','none', ...
	'Name','DP_TT Browser:  marquardt_torsten', ...
	'NumberTitle','off', ...
	'PaperType','a4letter', ...
	'PointerShapeCData',mat1, ...
	'Position',[0 0 476 439], ...
	'Resize','off', ...
	'Tag','Fig4', ...
	'UserData',mat2);
b = uicontrol('Parent',a, ...
	'Callback','dp_tt_browser new', ...
	'Position',[10 110 50 20], ...
	'String','New', ...
	'Tag','Pushbutton1');
b = uicontrol('Parent',a, ...
	'Callback','dp_tt_browser view', ...
	'Position',[10 85 50 20], ...
	'String','View', ...
	'Tag','Pushbutton1', ...
	'Value',1);
b = uicontrol('Parent',a, ...
	'Callback','dp_tt_browser delete', ...
	'Position',[10 60 50 20], ...
	'String','Delete', ...
	'Tag','Pushbutton2');
b = uicontrol('Parent',a, ...
	'Callback','dp_tt_browser probe_fit', ...
	'Position',[10 35 50 20], ...
	'String','probe fit', ...
	'Tag','Pushbutton2');
b = uicontrol('Parent',a, ...
	'Callback','dp_tt_browser help', ...
	'Position',[10 5 50 20], ...
	'String','Help', ...
	'Tag','Pushbutton2');
b = uicontrol('Parent',a, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','left', ...
	'Position',[7 191 55 15], ...
	'String','Created by:', ...
	'Style','text', ...
	'Tag','StaticText1');
b = uicontrol('Parent',a, ...
	'BackgroundColor',[1 1 1], ...
	'Position',[65 191 115 15], ...
	'Style','text', ...
	'Tag','created');
b = uicontrol('Parent',a, ...
	'BackgroundColor',[1 1 1], ...
	'Position',[253 191 27 15], ...
	'Style','text', ...
	'Tag','t_avg');
b = uicontrol('Parent',a, ...
	'BackgroundColor',[1 1 1], ...
	'Callback','dp_tt_browser comment', ...
	'HorizontalAlignment','left', ...
	'Max',12, ...
	'Position',[64 6 410 181], ...
	'Style','edit', ...
	'Tag','comment');
b = uicontrol('Parent',a, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','left', ...
	'Position',[198 192 51 15], ...
	'String','Averages:', ...
	'Style','text', ...
	'Tag','StaticText3');
b = uicontrol('Parent',a, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','left', ...
	'Position',[12 168 49 20], ...
	'String','comment:', ...
	'Style','text', ...
	'Tag','StaticText1');
b = uicontrol('Parent',a, ...
	'BackgroundColor',[1 1 1], ...
	'Callback','dp_tt_browser fill_text', ...
	'FontName','Courier', ...
	'Position',[5 211 472 213], ...
	'Style','listbox', ...
	'Tag','series', ...
	'Value',1);
b = uicontrol('Parent',a, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','right', ...
	'Position',[175 424 79 13], ...
	'String','Default Name:', ...
	'Style','text', ...
	'Tag','StaticText2');
b = uicontrol('Parent',a, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','left', ...
	'Position',[8 425 79 12], ...
	'String','Date:', ...
	'Style','text', ...
	'Tag','StaticText2');
b = uicontrol('Parent',a, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','right', ...
	'Position',[430 191 9 15], ...
	'String','#', ...
	'Style','text', ...
	'Tag','StaticText1');
b = uicontrol('Parent',a, ...
	'BackgroundColor',[1 1 1], ...
	'Position',[446 191 26 15], ...
	'Style','text', ...
	'Tag','n_meas');
b = uicontrol('Parent',a, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','left', ...
	'Position',[285 188 21 17], ...
	'String','[s]', ...
	'Style','text', ...
	'Tag','StaticText3');
b = uicontrol('Parent',a, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','left', ...
	'Position',[336 189 21 17], ...
	'String','Ear:', ...
	'Style','text', ...
	'Tag','StaticText3');
b = uicontrol('Parent',a, ...
	'BackgroundColor',[1 1 1], ...
	'Position',[359 191 38 15], ...
	'Style','text', ...
	'Tag','ear');

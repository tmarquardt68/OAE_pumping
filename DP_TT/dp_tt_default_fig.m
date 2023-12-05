function fig = dp_tt_default_fig()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
% This problem is solved by saving the output as a FIG-file.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.
% 
% NOTE: certain newer features in MATLAB may not have been saved in this
% M-file due to limitations of this format, which has been superseded by
% FIG-files.  Figures which have been annotated using the plot editor tools
% are incompatible with the M-file/MAT-file format, and should be saved as
% FIG-files.

load dp_tt_default_fig

h0 = figure('Color',[0.75 0.75 0.75], ...
	'Colormap',mat0, ...
	'MenuBar','none', ...
	'Name','DP_TT Default Browser', ...
	'NumberTitle','off', ...
	'PaperPosition',[18 180 576 432], ...
	'PaperType','a4letter', ...
	'PaperUnits','points', ...
	'Position',[0 3 526 395], ...
	'Resize','off', ...
	'Tag','Fig4');
h1 = uicontrol('Parent',h0, ...
	'Callback','dp_tt_default load', ...
	'ListboxTop',0, ...
	'Position',[23 8 53 20], ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ...
	'String','Load', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'Callback','dp_tt_default delete', ...
	'ListboxTop',0, ...
	'Position',[23 31 53 20], ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ... 
	'String','Delete', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[290 266 57 19], ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ... 
	'String','Averages:', ...
	'Style','text', ...
	'Tag','StaticText2');
h1 = uicontrol('Parent',h0, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[3 195 100 20], ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ... 
	'String','Frequencies f2 [Hz]:', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[3 171 100 20], ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ... 
	'String','Levels l2 (<75dB):', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[3 147 100 20], ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ... 
	'String','Freq. f1 = f2 / x ;  x=:', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[3 123 100 20], ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ... 
	'String','Lev. l1 = l2 + x ;  x=:', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[351 271 33 17], ...
	'Style','text', ...
	'Tag','n_avg');
h1 = uicontrol('Parent',h0, ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[108 174 409 17], ...
	'Style','text', ...
	'Tag','L2');
h1 = uicontrol('Parent',h0, ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[108 150 409 17], ...
	'Style','text', ...
	'Tag','F1');
h1 = uicontrol('Parent',h0, ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[108 126 409 17], ...
	'Style','text', ...
	'Tag','L1');
h1 = uicontrol('Parent',h0, ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[108 268 41 18], ...
	'Style','text', ...
	'Tag','created');
h1 = uicontrol('Parent',h0, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[4 265 100 20], ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ...
	'String','Created by:', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[108 246 409 17], ...
	'Style','text', ...
	'Tag','F_TT');
h1 = uicontrol('Parent',h0, ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[108 222 409 17], ...
	'Style','text', ...
	'Tag','L_TT');
h1 = uicontrol('Parent',h0, ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[108 198 409 17], ...
	'Style','text', ...
	'Tag','F2');
h1 = uicontrol('Parent',h0, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[2 219 100 20], ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ... 
	'String','Levels TT(<115dB):', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[2 243 100 20], ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ... 
	'String','Frequencies TT [Hz]:', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ...
	'BackgroundColor',[1 1 1], ...
	'Callback','dp_tt_default fill', ...
	'Position',[106 317 413 76], ...
	'Style','listbox', ...
	'Tag','defaults', ...
	'Value',1);
h1 = uicontrol('Parent',h0, ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ...
	'BackgroundColor',[1 1 1], ...
	'Callback','dp_tt_default comment', ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Max',12, ...
	'Position',[107 6 411 69], ...
	'Style','edit', ...
	'Tag','purpose');
h1 = uicontrol('Parent',h0, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[52 55 49 20], ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ... 
	'String','Comment:', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[416 266 75 20], ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ... 
	'String','Nested:', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[496 271 18 18], ...
	'Style','text', ...
	'Tag','nested_mode');
h1 = uicontrol('Parent',h0, ...
	'FontUnits','pixels', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','dp_tt_default rename', ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[107 293 411 20], ...
	'Style','edit', ...
	'Tag','rename');
h1 = uicontrol('Parent',h0, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[5 288 98 20], ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ... 
	'String','rename:', ...
	'Style','text', ...
	'Tag','StaticText1');

h1 = uicontrol('Parent',h0, ...
	'BackgroundColor',[0.752941 0.752941 0.752941], ...
	'Units','pixels', ...
	'ListboxTop',0, ...
	'Position',[30 360 75 20], ...
	'FontUnits','pixels', ... 
	'FontSize',10.6667, ... 
	'String','Default list:', ...
	'Style','text', ...
	'Tag','StaticText3');

if nargout > 0, fig = h0; end
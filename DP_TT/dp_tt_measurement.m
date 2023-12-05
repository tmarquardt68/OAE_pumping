function DP_TT_measurement(action, patient, series)
%     name:   Is either the default file name (if $action = 'load_defaults')
%             or the patient name (if $action = 'ini')
%
% Part of the OAE toolbox
% Copyright (C) 2008 Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).
%
% MODIFICATION:
% - writes into file, reads the file averaged and saves it
%   the old way (but named $avgN).
% - matrix $atomic is replaced by several $avgN variables
%   whre N is a 1,2,3... (atomic as a matrix became to big in cases of
%   long measurement series.)

global SAMPLE_RATE MIN_FREQ OAE_PATH

switch action
    case 'ini'
        h = findobj('Tag','dp_tt_measurement_fig');
        if isempty(h)
            h = dp_tt_measurement_fig;
        else
            figure(h)
        end
        try
            load([OAE_PATH,'\Subjects\' patient '/' series],...
                'F1_str','F2_str','F_TT_str','L1_str','L2_str','L_TT_str','nested_mode',...
                't_avg','comment')
            set(findobj(h,'Tag','F1'),'String',F1_str),
            set(findobj(h,'Tag','F2'),'String',F2_str),
            set(findobj(h,'Tag','F_TT'),'String',F_TT_str),
            set(findobj(h,'Tag','L1'),'String',L1_str),
            set(findobj(h,'Tag','L2'),'String',L2_str),
            set(findobj(h,'Tag','L_TT'),'String',L_TT_str),
            set(findobj(h,'Tag','nested_mode'),'Value',nested_mode),
            set(findobj(h,'Tag','t_avg'),'String',t_avg),
            set(findobj(h,'Tag','comment'),'String',comment),
        catch
        end
        set(gcf,'Name',['DP_TT Measurement: ',patient]),

    case 'load_defaults'
        % fill boxes with defaults
        h = findobj('Tag','dp_tt_measurement_fig');
        if exist([OAE_PATH,'\DP_TT\DP_TT_default_data\',patient],'file')
            load([OAE_PATH,'\DP_TT\DP_TT_default_data\',patient])
        else
            F1='';F2='';L1='';L2='';L_TT='';F_TT='';nested_mode=0;
        end
        try
            set(findobj(h,'Tag','F1'),'String',F1),
            set(findobj(h,'Tag','F2'),'String',F2),
            set(findobj(h,'Tag','F_TT'),'String',F_TT),
            set(findobj(h,'Tag','L1'),'String',L1),
            set(findobj(h,'Tag','L2'),'String',L2),
            set(findobj(h,'Tag','L_TT'),'String',L_TT),
            set(findobj(h,'Tag','nested_mode'),'Value',nested_mode),
            set(findobj(h,'Tag','t_avg'),'String',t_avg),
            set(findobj(h,'Tag','comment'),'String',comment),
            set(findobj(h,'Tag','defaults'),'String',default),
        catch
            set(findobj(h,'Tag','F1'),'String',F1_str),
            set(findobj(h,'Tag','F2'),'String',F2_str),
            set(findobj(h,'Tag','F_TT'),'String',F_TT_str),
            set(findobj(h,'Tag','L1'),'String',L1_str),
            set(findobj(h,'Tag','L2'),'String',L2_str),
            set(findobj(h,'Tag','L_TT'),'String',L_TT_str),
            set(findobj(h,'Tag','nested_mode'),'Value',nested_mode),
            set(findobj(h,'Tag','t_avg'),'String',t_avg),
            set(findobj(h,'Tag','comment'),'String',comment),
            set(findobj(h,'Tag','defaults'),'String',default),

        end


	case 'browse'
		dp_tt_default ini,

	case 'calibration'
		dp_tt_calibration('ini'),

	case 'measure'
		name = get(gcbf,'Name');
		if strcmp(get(gcbo,'String'),'Fit Tones: On'),
			set(gcbo,'String','Probe Fit: Off')
		end,
		measure(name(20:length(name))) % see function 'measure()' below!

	case 'close'


	case 'probe_fit',
		% plays from all three reciver to give audible feeback during probefit
        l = SAMPLE_RATE/MIN_FREQ; % longest possible period

		if strcmp(get(gcbo,'String'),'Fit Tones: On'),
			set(gcbo,'String','Probe Fit: Off'),
			drawnow,
            h_sound_device = sound_io('open_device',SAMPLE_RATE, 4, 0);
			while (1),
				sound_io('play',h_sound_device, ...
                    [sin(linspace(0,1000/MIN_FREQ*2*pi*(1-1/l), 3*l));zeros(3,3*l)]'/10);
                pause(.3)
				if strcmp(get(gcbo,'String'),'Fit Tones: On'),
					break,
				end,
				sound_io('play',h_sound_device, ...
                    [zeros(1,3*l);sin(linspace(0,1000/MIN_FREQ*2*pi*(1-1/l), 3*l));zeros(2,3*l)]'/10);
                pause(.3)
				if strcmp(get(gcbo,'String'),'Fit Tones: On'),
					break,
				end,
				sound_io('play',h_sound_device, ...
                    [zeros(2,3*l);sin(linspace(0,1000/MIN_FREQ*2*pi*(1-1/l), 3*l));zeros(1,3*l)]'/10);
				pause(.3)
				if strcmp(get(gcbo,'String'),'Fit Tones: On'),
					break,
				end,
			end,
            sound_io('close_device',h_sound_device);
		else
			set(gcbo,'String','Fit Tones: On')
		end,


	case 'ear_l',
		set(findobj(gcbf,'Tag','ear_r'),'Value',0),

	case 'ear_r',
		set(findobj(gcbf,'Tag','ear_l'),'Value',0),

	case 'help'
		return

	otherwise
		error('Unkown action (case)!');
end %switch action

%======================================================================

function measure(patient)
% calls dp_tt_probe_check.m


global SAMPLE_RATE MIN_FREQ OAE_PATH BITS_PER_SAMPLE

global CURRENT_EAR % filled by dp_tt_probe_check.m
%  CURRENT_EAR(:,1) - rec1 transfer fct. (complex)
%  CURRENT_EAR(:,2) - rec2 transfer fct. (complex)
%  CURRENT_EAR(:,3) - TT   transfer fct. (complex)
%  CURRENT_EAR(:,4) - mic transfer fct. (complex)
%  CURRENT_EAR(1,5) - max. soundpressure level of rec1 [dB SPL]
%  CURRENT_EAR(2,5) - max. soundpressure level of rec2 [dB SPL]
%  CURRENT_EAR(3,5) - max. soundpressure level of TT [dB SPL]
%  CURRENT_EAR(:,6) - transfer fct.(complex) of current calibration


load([OAE_PATH ,...
    '\DP_TT\DP_TT_calibration_data\00000_last_calibration'],'electrostat')
l = SAMPLE_RATE/MIN_FREQ;

set(findobj('Tag','probe_fit'),'String','Fit Tones: On');
drawnow,


% get prarameter from figure
h_fig = gcbf;
F1_str = get(findobj(h_fig,'Tag','F1'),'String');
F2_str = get(findobj(h_fig,'Tag','F2'),'String');
F_TT_str = get(findobj(h_fig,'Tag','F_TT'),'String');
L1_str = get(findobj(h_fig,'Tag','L1'),'String');
L2_str = get(findobj(h_fig,'Tag','L2'),'String');
L_TT_str = get(findobj(h_fig,'Tag','L_TT'),'String');
t_avg = eval(get(findobj(h_fig,'Tag','t_avg'),'String'));
comment = get(findobj(h_fig,'Tag','comment'),'String');
default = get(findobj(h_fig,'Tag','defaults'),'String');
nested_mode = get(findobj(h_fig,'Tag','nested_mode'),'value');
created = get(findobj(h_fig,'Tag','created'),'String');
ear_l = get(findobj(h_fig,'Tag','ear_l'),'Value');
ear_r = get(findobj(h_fig,'Tag','ear_r'),'Value');
probe_check = get(findobj(h_fig,'Tag','probe_check'),'Value');


% save as default
if ~exist([OAE_PATH,'\DP_TT\DP_TT_default_data']),
	eval(['!mkdir ', OAE_PATH,'\DP_TT\DP_TT_default_data'])
end
sample_rate = SAMPLE_RATE;
min_freq = MIN_FREQ;
save([OAE_PATH ,'\DP_TT\DP_TT_default_data\',default],...
	'F1_str','F2_str','F_TT_str','L1_str','L2_str','L_TT_str', 'nested_mode',...
	't_avg','created','comment','default','min_freq'),
save([OAE_PATH ,'\DP_TT\DP_TT_default_data\_last_default'],...
	'F1_str','F2_str','F_TT_str','L1_str','L2_str','L_TT_str', 'nested_mode',...
	't_avg','created','comment','default','min_freq'),

L1 = eval(['[',L1_str,']']);
L2 = eval(['[',L2_str,']']);
L_TT = eval(['[',L_TT_str,']']);
F1 = eval(['[',F1_str,']']);
F2 = eval(['[',F2_str,']']);
if F_TT_str(1) == 'n' | F_TT_str(1) == 'N'
    F_TT_str = inf;
else
    F_TT = eval(['[',F_TT_str,']']);
    if min(F_TT)< MIN_FREQ,
        msgbox(['Frequency TT must be at least ' num2str(2*MIN_FREQ) ' Hz !'],...
            'ERROR - DP_TT Measurement:','error');
        return,
    end
%     if sum(mod(F_TT,MIN_FREQ))>0,
%         msgbox(strvcat('At least one of Frequencies TT has been rounded',...
%             ['to a multiple of ' num2str(MIN_FREQ) ' Hz !']),...
%             'WARNING - DP_TT Measurement:','warn');
%     end
end
% some checks of the parameter
% if length(F1)==0|length(L1)==0|length(F2)==0|length(L2)==0|...
% 		length(F_TT)==0|length(L_TT)==0|length(t_avg)==0|...
% 		length(created)==0|(ear_l==0 & ear_r==0)|length(L_TT)==0,
    if length(F1)==0|length(L1)==0|length(F2)==0|length(L2)==0|...
		length(F_TT)==0|length(L_TT)==0|length(t_avg)==0|...
		length(L_TT)==0
	msgbox('You left some parameter for the measurement undefined!',...
		'ERROR - DP_TT Measurement:','error');
	return,
end

% if sum(mod(F1,MIN_FREQ))>0,
% 	msgbox(strvcat('At least one of Frequencies f1 has been rounded',...
% 		['to a multiple of ' num2str(MIN_FREQ) ' Hz !']),...
% 		'WARNING - DP_TT Measurement:','warn');
% end
% if sum(mod(F2,MIN_FREQ))>0,
% 	msgbox(strvcat('At least one of Frequencies f2 has been rounded',...
% 		['to a multiple of ' num2str(MIN_FREQ) ' Hz !']),...
% 		'WARNING - DP_TT Measurement:','warn');
% end


n = 0;n1=0;
if nested_mode,
    F_1=F1; F_2=F2; F__TT=F_TT; L_1=L1; L_2=L2; L__TT=L_TT;
    for f1=F_1
        for f2=F_2
            for f_TT=F__TT
                for l2=L_2
                    for l1=L_1
                        for l_TT=L__TT
                            n = n+1;
                            F1(n)=f1;
                            F2(n)=f2;
                            F_TT(n)=f_TT;
                            L1(n)=l1;
                            L2(n)=l2;
                            L_TT(n)=l_TT;
                        end
                    end
                end
            end
        end
    end
end

% number of single measurements:
n=max([length(F1),length(F2),length(F_TT),length(L1),length(L2),...
	length(L_TT)]);
% fill parameter up to the same length
if length(F1)==1,F1 = F1*ones(1,n); end,
if length(F2)==1,F2 = F2*ones(1,n); end,
if length(F_TT)==1,F_TT = F_TT*ones(1,n); end,
if length(L1)==1,L1 = L1*ones(1,n); end,
if length(L2)==1,L2 = L2*ones(1,n); end,
if length(L_TT)==1,L_TT = L_TT*ones(1,n); end,
% convert nested to non-nested parameter

if length(F1)~=n|length(L1)~=n|length(F2)~=n|length(L2)~=n...
		|length(F_TT)~=n|length(L_TT)~=n,
	msgbox(strvcat('ERROR - DP_TT_Measurement:',...
		'In non-nested mode all levels and frequencies ',...
		'must be of same length or scalar!'));
	return,
end

% convert freqenciies from Hz to spectral line number
F2 = round(F2);
F1 = round(F2./F1);
L1 = L2 + L1;
L1 = round(L1);
L2 = round(L2);

% convert freqenciies from Hz to spectral line number
F2 = 2*round(F2/2/MIN_FREQ)-1;
F1 = 2*round(F1/2/MIN_FREQ)-1;
% F_TT = 2*round(F_TT/2/MIN_FREQ);
F_TT = F_TT/MIN_FREQ;

if probe_check == 0, % CURRENT_EAR is filled with data of last calibration
    % show measurement time
    if ~strcmp('Yes',questdlg(strvcat(strvcat(...
            'Time for this Measurement series will be about: ',...
            [num2str(floor(n*t_avg/60)),' min  ', ...
            num2str(round(mod(n*t_avg,60))),' sec  '...
            num2str((F2(1))*MIN_FREQ) ' ' num2str((F1(1))*MIN_FREQ) ' ' ...
            num2str((2*F1(1)-F2(1))*MIN_FREQ) ' ' num2str((F2(1)-F1(1))*MIN_FREQ) ' ' ...
            num2str((F2(1)-F1(1)./F_TT(1))*MIN_FREQ)],...
            'Press "Yes" to continue!')),'MESSAGE: DP_TT measerment:'))
        return,
    end
    % load mat-file containing last calibration
    dp_tt_probe_check('load_last');
else
    % show measurement time
    if ~strcmp('Yes',questdlg(strvcat(strvcat(...
            'Time for this Measurement series will be about: ',...
            [num2str(floor(n*t_avg/60)),' min  ', ...
            num2str(round(mod(n*t_avg,60))),' sec  + 30 sec   '...
            num2str((F2(1))*MIN_FREQ) ' ' num2str((F1(1))*MIN_FREQ) ' ' ...
            num2str((2*F1(1)-F2(1))*MIN_FREQ) ' ' num2str((F2(1)-F1(1))*MIN_FREQ) ' ' ...
            num2str((F2(1)-F1(1)./F_TT(1))*MIN_FREQ)],...
            'Press "Yes" to continue!')),'MESSAGE: DP_TT measerment:'))
        return,
    end
    dp_tt_probe_check('ini'),
end

if isempty(CURRENT_EAR),
	msgbox(strvcat('ERROR - DP_TT Measurement:',...
		'Probe check failed!'));
	return,
end,
pause(0.5),


Hdb_rec1 = 20*log10(abs(CURRENT_EAR(F1+1,1)'));
Hdb_rec2 = 20*log10(abs(CURRENT_EAR(F2+1,2)'));
Hdb_TT = 20*log10(abs(CURRENT_EAR(F_TT+1,3)'));
if isinf(F_TT(1))
    Hdb_TT = inf;
else
    Hdb_TT = 20*log10(abs(CURRENT_EAR(F_TT+1,3)'));
    if max(L_TT - Hdb_TT) > CURRENT_EAR(3,5)
        msgbox('The corrected Level L_TT is higher than max_level!',...
            'ERROR - DP_TT Measurement:','error');
        return,
    end
end
if max(L2 - Hdb_rec2) > CURRENT_EAR(2,5)
	msgbox('The corrected Level L2 is higher than max_level!',...
		'ERROR - DP_TT Measurement:','error');
	return,
end
if max(L1 - Hdb_rec1) > CURRENT_EAR(1,5)
	msgbox('The corrected Level L1 is higher than max_level!',...
		'ERROR - DP_TT Measurement:','error');
	return,
end


% open new series file
series = ['dp_tt',datestr(now,30)];
if ~exist([OAE_PATH,'\Subjects\',patient],'file')
	eval(['!mkdir ', OAE_PATH,'\Subjects\',patient])
end

% save some data of series
% H_TT/H_rec1 of current ear / calibration (normalized) = probe_fit
HdB_ear_channel = 20*log10(...
    CURRENT_EAR(1:round(10000/MIN_FREQ)+1,1))*...
	abs(CURRENT_EAR(round(1000/MIN_FREQ)+1,6))./...
	CURRENT_EAR(1:round(10000/MIN_FREQ)+1,6);
[m, b] = size(comment);
comment = [[default,zeros(1,b-length(default))];...
	[comment,zeros(m,length(default)-b)]];
l1 = []; l2 = []; l_TT = []; l_dp = [];
h_mic = CURRENT_EAR(:,4);
save([OAE_PATH,'\Subjects\',patient,'\',series],...
	't_avg','created','comment','nested_mode','ear_l','patient',...
	'min_freq','sample_rate','HdB_ear_channel','h_mic',...
    'F1_str','F2_str','F_TT_str','L1_str','L2_str','L_TT_str',...
    'F1','L1','l1','F2','L2','l2','F_TT','L_TT','l_TT','l_dp'),

% open new viewer figure
dp_tt_viewer('ini',patient,series),

% start of measurement series
n_avg = t_avg*MIN_FREQ;
SecToRecord=(n_avg+2)*l/SAMPLE_RATE;
PsychPortAudio('Close')
h_sound_device = sound_io('open_device',SAMPLE_RATE, 3, 2);
PsychPortAudio('GetAudioData', h_sound_device, 2*SecToRecord);
PsychPortAudio('FillBuffer',h_sound_device,repmat(zeros(3,l),1,n_avg+2));
PsychPortAudio('Start', h_sound_device);
l1=NaN*ones(1,n);l2=l1;l_TT=l1;timeStamps=l1;
a = 1;
t_loop = tic;
while a <= n+1
    if a<=n
    	% phase correction
    	wave = zeros(3,l);
    	[re, im] = pol2cart(-pi/2-angle(CURRENT_EAR(F1(a)+1,1)),l/2);
    	wave(1,F1(a)+1) = re + im*1i;
    	wave(1,l - F1(a)+1) = re - im*1i;
    	[re, im] = pol2cart(-pi/2-angle(CURRENT_EAR(F2(a)+1,2)),l/2);
    	wave(2,F2(a)+1) = re + im*1i;
    	wave(2,l - F2(a)+1) = re - im*1i;
        wave(1,:) = real(ifft(wave(1,:)));
        wave(2,:) = real(ifft(wave(2,:)));


        ogain1 = 10^( (L1(a)-Hdb_rec1(a)-CURRENT_EAR(1,5))/20);
        ogain2 = 10^( (L2(a)-Hdb_rec2(a)-CURRENT_EAR(2,5))/20);
        wave_out(1,:)=wave(1,:)*ogain1;
        wave_out(2,:)=wave(2,:)*ogain2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if electrostat % for electrostatic speakers only
            wave_out = ((200*wave_out+212.7).^0.5-sqrt(212.7))/sqrt(200);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if L_TT(a) > -90
            [re, im] = pol2cart(-pi/2-angle(CURRENT_EAR(F_TT(a)+1,1)),l/2);
            [re, im] = pol2cart(-pi/2,l/2);
            wave(3,F_TT(a)+1) = re + im*1i;
            wave(3,l - F_TT(a)+1) = re - im*1i;
            wave(3,:) = real(ifft(wave(3,:)));
            ogainT = 10^( (L_TT(a)-Hdb_TT(a)-CURRENT_EAR(3,5))/20);
            wave_out(3,:) = wave(3,:)*ogainT;
        else
            wave_out(3,:) = zeros(size(wave(3,:)));
        end
        % Volume setting via wave amplitude!
    end
	% the measurement
    t_rec=tic;
    recBuffer = PsychPortAudio('GetAudioData',h_sound_device, [],SecToRecord,SecToRecord)';
    recording_time =toc(t_rec);
    pause(0.1)
    loop_time = toc(t_loop)
    if loop_time > 6 
        h = findobj('Tag','dp_tt_viewer_fig');
        set(h(end),'Color','r')
    end
    PsychPortAudio('RefillBuffer',h_sound_device,0,repmat(wave_out,1,n_avg+2));
    PsychPortAudio('Start', h_sound_device);
    t_loop = tic;

    if a>1
        a= a-1;
        %average
        avg_wave_in = double(squeeze(sum(reshape(recBuffer(1*l+1:(n_avg+1)*l,:), l, n_avg,[]), 2)./n_avg));

    	H_avg_acoust = fft(avg_wave_in(:,1))./h_mic;
    	H_avg_CM = fft(avg_wave_in(:,2));

        % get the actual levels
        f(1) = 2*F1(a)-F2(a)+1;
        f(2) = F2(a)-F1(a)+1;
        f(3) = 2*F1(a)-F2(a)-2*F_TT(a)+1;
        f(4) = 2*F1(a)-F2(a)-F_TT(a)+1;
        f(5) = 2*F1(a)-F2(a)+F_TT(a)+1;
        f(6) = 2*F1(a)-F2(a)+2*F_TT(a)+1;
        f(7) = F2(a)-F1(a)-2*F_TT(a)+1;
        f(8) = F2(a)-F1(a)-F_TT(a)+1;
        f(9) = F2(a)-F1(a)+F_TT(a)+1;
        f(10) = F2(a)-F1(a)+2*F_TT(a)+1;
        for q=length(f):-1:1
            if f(q)>0 & f(q) < length (H_avg_acoust)
                l_dp(a,q) = round(20*log10(abs(H_avg_acoust(f(q))))*10)/10;
                l_dp(a,q+10) = round(20*log10(abs(H_avg_CM(f(q))))*10)/10;
                phase_dp(a,q) = unwrap(angle(H_avg_acoust(f(q))))/2/pi;
                phase_dp(a,q+10) = unwrap(angle(H_avg_CM(f(q))))/2/pi;
            else
                l_dp(a,q) = nan;
                l_dp(a,q+10) = nan;
                phase_dp(a,q) = nan;
                phase_dp(a,q+10) = nan;
            end
        end
        l_TT(a) = round(20*log10(abs(H_avg_acoust(F_TT(a)+1)))*10)/10;
        l1(a) = round(20*log10(abs(H_avg_acoust(F1(a)+1))));
        l2(a) = round(20*log10(abs(H_avg_acoust(F2(a)+1))));


    	% save data
    	eval(['avg', num2str(a), '= avg_wave_in;']);
        timeStamps(a)= now*(24*60*60); % in seconds
        save([OAE_PATH,'\Subjects\',patient,'\',series],...
            ['avg', num2str(a)],'timeStamps','a','F1','F2','F_TT','L1','L2','L_TT',...
            'l1','l2','l_TT','l_dp','phase_dp','-append'),

        % show results
        dp_tt_viewer('fill',patient,series), % update viewer figure
        dp_tt_viewer('plot_replace',patient,series),
        a=a+1;
    end

% 	% get new L_TT value from user,
% 	new_L_TT = 0;
% 	while ~isempty(new_L_TT) & new_L_TT == 0,
% 		new_L_TT = inputdlg(strvcat(...
% 			['Give new L_TT level (old was: ',num2str(L_TT(a)),' dB).'], ...
% 			'Press OK (leave the same value) to repeat measurment.', ...
% 			'Press CANCEL (leave blank) if last is ok.',...
% 			'Give (different) negative level to terminate meauerment.'),...
% 			'DP_TT Measurement', 1, {num2str(L_TT(a))});
% 
%         new_L_TT = '';
% 
% 		if ~isempty(new_L_TT),
% 
% 			try,
% 				new_L_TT = str2num(new_L_TT{1});
% 			catch,
% 				new_L_TT = 0;
% 			end,
% 
% 			if length(new_L_TT)~=1,
% 				uiwait(msgbox(strvcat('Please enter a number',...
% 					'Use a decimal POINT !'),...
% 					'ERROR DP_TT_measurement:','error'));
% 				new_L_TT = 0;
% 			end,
% 
% 			if max(new_L_TT - Hdb_TT(a)) > CURRENT_EAR(1,5),
% 				h = msgbox('The corrected Level L_TT is higher than max_level!',...
% 					'ERROR - DP_TT Measurement:','error');
% 				uiwait(h),
% 				new_L_TT = 0;
% 			end,
% 		end, % ~isempty(new_L_TT)
% 	end, % while(new_L_TT == 0)
% 	drawnow,
% 	
% 	if isempty(new_L_TT), % do next parameter set
% 		clear(['avg', num2str(a)]),
% 		a = a+1;
% 	else, % insert new measurement with new_L_TT
% 		if new_L_TT < -1; % stop measurement
% 			return,
%         elseif new_L_TT == -1
%             break
% 		end,
% 		F1 = [F1(1:a), F1(a:n)];
% 		L1 = [L1(1:a), L1(a:n)];
% 		F2 = [F2(1:a), F2(a:n)];
% 		L2 = [L2(1:a), L2(a:n)];
% 		F_TT = [F_TT(1:a), F_TT(a:n)];
% 		L_TT = [L_TT(1:a), new_L_TT, L_TT(a+1:n)+new_L_TT-L_TT(a)];
% 		Hdb_rec1 = [Hdb_rec1(1:a), Hdb_rec1(a:n)];
% 		Hdb_rec2 = [Hdb_rec2(1:a), Hdb_rec2(a:n)];
% 		Hdb_TT = [Hdb_TT(1:a), Hdb_TT(a:n)];
% 		n = max(length(F1));
% 		clear(['avg', num2str(a)]),
% 		a = a+1;
%     end,    

    clear(['avg', num2str(a)]),
    a = a+1;
    dp_tt_browser('fill',patient),
end, % while(a <= n),
sound_io('close_device',h_sound_device);
close(gcbf);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


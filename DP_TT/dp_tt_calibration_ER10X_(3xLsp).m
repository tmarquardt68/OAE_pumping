% DP_TT_CALIBRATION(action):
% Does the calibration and writes calibration data in the directory 
% ...\OAE\DP_TT\DP_TT_calibration_data as a mat-file.
%
% Part of the OAE toolbox 
% Copyright (C) 2008 Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).

% Saves the following variables in a mat-file 
% (e.g. '00000_last_calibration.mat' used during the 
%  following DP_TT measurements)
%
%  H_rec1 - rec1 transfer fct. (complex)
%  H_rec2 - rec2 transfer fct. (complex)
%  H_TT   - TT transfer fct. (complex)
%  H_mic  - mic transfer fct. (complex)
%  max_level_rec1 - max. soundpressure level of rec1 [dB SPL]
%  max_level_rec2 - max. soundpressure level of rec2 [dB SPL]
%  max_level_rec3) - max. soundpressure level of TT [dB SPL]
%  mic_level_correction - mic correction level at 1kHz [dB]
%  snr_rec1 - SNR + THD test spectrum of rec1  [dB SPL]
%  snr_rec2 - SNR + THD test spectrum of rec2  [dB SPL]
%  snr_TT   - SNR + THD test spectrum of TT    [dB SPL]

function dp_tt_calibration(action)

global MIN_FREQ SAMPLE_RATE OAE_PATH

l = SAMPLE_RATE/MIN_FREQ; % longest possible period
TT_ref_freq = 10;

if ~exist('action'),
    action = 'ini';
end,

switch action
    case 'ini',
        if (~exist([OAE_PATH,'\DP_TT\DP_TT_calibration_data']))
			eval(['!mkdir ', OAE_PATH,'\DP_TT\DP_TT_calibration_data'])
		end
% 		try,
% 			delete([OAE_PATH,...
% 				'\DP_TT\DP_TT_calibration_data\00000__NONE__.mat']),
% 		end,
		h_calib = findobj('Tag','dp_tt_calibration_fig');
		if(isempty(h_calib)),
			h_calib = dp_tt_calibration_fig;
		else,
			figure(h_calib)
        end,
        h_ref = findobj(gcf,'Tag','references');
        set(h_ref,'Value',1);
		setappdata(h_calib,'BK_amplitude_94dBspl', []), % forces 
        dp_tt_calibration('fill'),
        
    case 'fill'    % fill listbox with references
        h_references = findobj(gcf,'Tag','references');
        str = getfield(what([OAE_PATH,...
                '\DP_TT\DP_TT_calibration_data']),'mat');
        if (isempty(str)),
            return;
        end,
        files(length(str),30)=' ';
        for(n=1:length(str)),
            files(n,1:length(str{n})) = str{n};
        end,
        str = cellstr(sortrows(files));
        set(h_references,'String',str),
        
        % fill reference comment box
        val = get(h_references,'Value');
        path_str = [OAE_PATH,...
                '\DP_TT\DP_TT_calibration_data\',str{val}];
        try,
            load (path_str,'comment'),
        catch,
            h = errordlg(['Unable to open file: ',path_str],...
                'ERROR - Calibration: data unit listbox');
            return,
        end,
        if (~exist('comment')),
            comment = {''};
        end,
        set(findobj(gcf,'Tag','ref_comment'),'String',comment),
        if strcmp(str{val},'00000__NONE__.mat')
            set(findobj(gcf,'Tag','comment'),'String',comment),
		end,
       
		
    case 'retry'
		t_avg = 20;
    	h_fig = findobj('Tag','dp_tt_calibration_fig');
		 
		% save data of previous attempt as last trial
        if exist([OAE_PATH,'\DP_TT\DP_TT_calibration_data\00000__NONE__.mat'])
            load ([OAE_PATH,...
				'\DP_TT\DP_TT_calibration_data\00000__NONE__.mat']),
		elseif exist([OAE_PATH,'\DP_TT\DP_TT_calibration_data\00000_last_trial.mat'])
            load ([OAE_PATH,...
				'\DP_TT\DP_TT_calibration_data\00000_last_trial.mat']),,
		else,
			tmp = zeros(1,l);
			H_rec1=tmp;H_rec2=tmp;H_TT=tmp;H_mic=tmp;snr_rec1=tmp;...
				snr_rec2=tmp; snr_TT=tmp;mic_level_correction=0;max_level_rec1=0;...
				comment=[];max_level_rec2=0;max_level_TT=0;min_freq=MIN_FREQ;...
				sample_rate=SAMPLE_RATE; BK_amplitude_94dBspl = -99;
		end,
		save([OAE_PATH,...
			'\DP_TT\DP_TT_calibration_data\00000_last_trial.mat'],...
			'H_rec1','H_rec2','H_TT','H_mic','snr_rec1','snr_rec2', ...
			'snr_TT','mic_level_correction','max_level_rec1','comment',...
			'max_level_rec2','max_level_TT','min_freq','sample_rate',...
			'BK_amplitude_94dBspl'),

		% get "numerical signal level" of BK mic when exposed to 94 dB SPL
		BK_amplitude_94dBspl = getappdata(h_fig,'BK_amplitude_94dBspl');
		if isempty(BK_amplitude_94dBspl),
			uiwait(msgbox(strvcat( ...
				'Connect the BK microphone with input channel 2 of AD converter',...
				'BK2636 :Set the gain on BK pre-amplifier to 20dB!',...
				'Put the BK microphone into the calibrator and switch it on'), ...
				'BK calibration')),

			h_sound_device = sound_io('open_device',SAMPLE_RATE, 4, 2);
            waveIn = sound_io('playrecord',h_sound_device,zeros(SAMPLE_RATE,4),1);
            % waveIn = sound_io('record', 3, 1); % record from channel 3 for x sec.
            sound_io('close_device',h_sound_device);



			BK_amplitude_94dBspl = mean(abs(waveIn(:,2)'));
			if ~strcmp('Yes',questdlg(strvcat( ...
					['BK_amplitude_94dBspl = ', num2str(BK_amplitude_94dBspl)],...
					'(should be around 0.5)'), ...
					'BK calibration')),
				return,
			end,
			setappdata(h_fig,'BK_amplitude_94dBspl',BK_amplitude_94dBspl);
			uiwait(msgbox(strvcat( ...
				'Connect BK microphone and OAE probe usinng 1 ccm test cavity.',...
				'Leave all setting of BK equipment.'), ...
				'DP_TT calibration')),
		end,
    
        set(h_fig,'Name','DP_TT Calibration in progress')
		drawnow
			
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SNR & THD test, probe input on AD channel 2
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tone = sin(linspace(0,1000/MIN_FREQ*2*pi*(1-1/l), l));
        snr_rec1 = data_acquisition(tone', 1, 2, ceil(t_avg*MIN_FREQ));
		snr_bk = snr_rec1(:,2);
		snr_rec1 = snr_rec1(:,1);
        snr_rec1 = fft(snr_rec1);
        
        snr_rec2 = data_acquisition([zeros(length(tone),1), tone'], 2, 1, ceil(t_avg*MIN_FREQ));
        snr_rec2 = fft(snr_rec2);
        
        tone = sin(linspace(0,1000/MIN_FREQ*2*pi*(1-1/l), l));
        snr_TT = data_acquisition([zeros(length(tone),2), tone'], 3, 1,ceil(t_avg*MIN_FREQ));      
        snr_TT = fft(snr_TT);
        
		% microphon correction
		amplitude_max_ouput = mean(abs(snr_bk));
		max_level_rec1 = 94 - 20*log10(BK_amplitude_94dBspl ./ amplitude_max_ouput);
		% In BK4157 use 93.3 dB SPL
        % max_level_rec1 = 93.3 - 20*log10(BK_amplitude_94dBspl ./ amplitude_max_ouput);
		mic_level_correction = max_level_rec1 - 20*log10(abs(snr_rec1(...
			round(1000/MIN_FREQ) + 1)));

		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TF measurements        
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% TF of probe mic
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% white noise generation
		rand('seed',0),
		theta = 2*pi*rand(1,l/2-1);
		[re im] = pol2cart([0 theta 0 ...
			-theta(l/2-1:-1:1)],...
			ones(1,l));
		H_wave = (re +im*i)';
		H_wave(1)=0;
				
		H_wave_HF = H_wave;
% 		H_wave_HF([1:round(1000/MIN_FREQ), ...
% 			l-round(1000/MIN_FREQ)+2:l])=0;
		wave_HF = real(ifft(H_wave_HF));
		wave_HF = wave_HF./max(abs(wave_HF));

		H_wave_LF = H_wave;
		H_wave_LF(round(1000/MIN_FREQ)+1:l-...
			round(1000/MIN_FREQ)+1)=0;
		wave_LF = real(ifft(H_wave_LF));
		wave_LF = wave_LF./max(abs(wave_LF));
		
		pause(3), % wait after overloading BK equipment with LF tone
		wave_in = data_acquisition([wave_HF,  zeros(length(wave_HF),2)],...
			3,2, ceil(max(20,t_avg)*MIN_FREQ));

		%load([OAE_PATH,'\DP_TT\H_bk']), %corection of BK calibration path
		H_mic = fft(wave_in(:,1))./fft(wave_in(:,2));%.*H_bk';%probe mic / corected B&K
		H_mic = H_mic/abs(H_mic(round(1000/...
			MIN_FREQ)+1))/10^(mic_level_correction/20);

		
		% TF of receivers
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		wave = real(ifft(H_wave));
		wave = wave./max(abs(wave));
		H_wave = fft(wave);

		% rec1
		H_rec1 = data_acquisition(wave, 1, 1, ceil(t_avg*MIN_FREQ));
		H_rec1 = fft(H_rec1)./H_wave./H_mic;
        snr_rec1 = 20*log10(abs(snr_rec1./H_mic));
        max_level_rec1 = snr_rec1(round(1000/MIN_FREQ)+1);

        % rec2
        H_rec2 = data_acquisition([zeros(length(wave),1), wave], 2, 1, ceil(t_avg*MIN_FREQ));
        H_rec2= fft(H_rec2)./H_wave./H_mic;
        snr_rec2 = 20*log10(abs(snr_rec2./H_mic));
        max_level_rec2 = snr_rec2(round(1000/MIN_FREQ)+1);
      
        % TT
        H_TT = data_acquisition([zeros(length(wave),2), wave], 3, 1, ceil(t_avg*MIN_FREQ));
        H_TT = fft(H_TT)./fft(wave)./H_mic;           
         snr_TT = 20*log10(abs(snr_TT./H_mic));
        max_level_TT = snr_TT(round(1000/MIN_FREQ)+1);
        
        % save calibration data
        comment = get(findobj(h_fig,'Tag','comment'),'String');
		[m n] = size(comment);
		comment = comment(2:m,:); % delete the date
        comment = [[datestr(now),zeros(1,n-length(datestr(now)))];comment];
		set(findobj(h_fig,'Tag','comment'),'String', comment);
		
        sample_rate = SAMPLE_RATE;
        min_freq = MIN_FREQ;
        save([OAE_PATH,...
            '\DP_TT\DP_TT_calibration_data\00000__NONE__.mat'],...
            'H_rec1','H_rec2','H_TT','H_mic','snr_rec1','snr_rec2', ...
            'snr_TT','mic_level_correction','max_level_rec1','comment',...
            'max_level_rec2','max_level_TT','min_freq','sample_rate',...
			'BK_amplitude_94dBspl'),
        
        set(h_fig, 'Name', 'DP_TT Calibration finished - now plotting')
        dp_tt_calibration('fill'),
        dp_tt_calibration('plot'),
        set(findobj(gcbf,'Tag','retry'),'String','Retry'),
        set(h_fig,'Name','DP_TT Calibration done. OK or RETRY?')

        
        
	case 'plot'
		% load reference and current data
		h_fig = findobj('Tag','dp_tt_calibration_fig');
		h_references = findobj(h_fig,'Tag','references');
		val = get(h_references, 'Value');
		str = get(h_references, 'String');
		if (isempty(str)),
			reference = '00000__NONE__';
		else,
			reference = str{val};
		end,
		load ([OAE_PATH,'\DP_TT\DP_TT_calibration_data\',reference]),
		[m n] = size(H_rec1);
		if m == 1,
			H_rec1=H_rec1'; H_rec2=H_rec2'; H_TT=H_TT'; H_mic=H_mic';
			snr_rec1=snr_rec1'; snr_rec2=snr_rec2'; snr_TT=snr_TT';
		end,
		ref_H_rec1 = H_rec1; ref_H_rec2 = H_rec2; ref_H_TT = H_TT;
		ref_H_mic = H_mic;ref_snr_rec1 =snr_rec1;ref_snr_rec2=snr_rec2;
		ref_snr_TT=snr_TT;ref_max_level_rec1=max_level_rec1;
		ref_max_level_rec2=max_level_rec2;ref_max_level_TT=max_level_TT;
		ref_mic_level_correction = mic_level_correction;

		if(exist([OAE_PATH,...
				'\DP_TT\DP_TT_calibration_data\00000__NONE__.mat'],'file')),
			load ([OAE_PATH,...
				'\DP_TT\DP_TT_calibration_data\00000__NONE__']),
		end,
		
		up_lim = round(20000/MIN_FREQ);
		lo_lim = round(995/MIN_FREQ);
		
		
		% plot phase calibration (transfer functions)
		x = linspace(0,SAMPLE_RATE/2,l/2+1);
		h = findobj('Tag','phase_fig');
		if (isempty(h)),
			phase_fig,
		else,
			figure(h),
		end,

        if length(H_mic) ~= length(ref_H_mic)
            H_mic = ref_H_mic;
            H_rec1 = ref_H_rec1;
            H_rec2 = ref_H_rec2;
            H_TT = ref_H_TT;
            snr_rec1 = ref_snr_rec1;
            snr_rec2 = ref_snr_rec2;
            snr_TT = ref_snr_TT;
            up_lim = length(ref_H_mic)/2
            x = linspace(0,SAMPLE_RATE/2,length(ref_H_mic)/2+1);
        end

            
		semilogx(x(3:up_lim),unwrap(angle(ref_H_mic(3:up_lim)))','g:',...
			x(lo_lim+2:up_lim),...
			unwrap(angle(ref_H_rec1(lo_lim+2:up_lim)))'+.025*x(lo_lim+2:up_lim),'r:',...
			x(lo_lim+2:up_lim),...
			unwrap(angle(ref_H_rec2(lo_lim+2:up_lim)))'+.025*x(lo_lim+2:up_lim),'b:',...
			x(3:lo_lim+1),...
			unwrap(angle(ref_H_TT(3:lo_lim+1)))','k:',...
			x(3:up_lim),...
			unwrap(angle(H_mic(3:up_lim)))','g',...
			x(lo_lim+2:up_lim),...
			unwrap(angle(H_rec1(lo_lim+2:up_lim)))'+ .55*x(lo_lim+2:up_lim),'r',...
			x(lo_lim+2:up_lim),...
			unwrap(angle(H_rec2(lo_lim+2:up_lim)))'+ .55*x(lo_lim+2:up_lim),'b',...
			x(3:lo_lim+1),...
			unwrap(angle(H_TT(3:lo_lim+1)))'+.55*x(3:lo_lim+1),'k.-'),
		grid on
		text(15,0,strcat('green = mic, red = rec1, blue = rec2,',...
			' black=TT, ---- =current, ....=reference)')),
		zoom on,

		% plot magnitude calibration (transfer functions)
		h = findobj('Tag','magnitude_fig');
		if (isempty(h)),
			magnitude_fig,
		else,
			figure(h),
		end,
		semilogx(x(2:up_lim),20*log10(abs(ref_H_mic(2:up_lim)))....
			+ ref_mic_level_correction,'g:',...
			x(2:up_lim),20*log10(abs(ref_H_rec1(2:up_lim)...
			/abs(ref_H_rec1(round(1000/MIN_FREQ)+1)))),'r:',...
			x(2:up_lim),20*log10(abs(ref_H_rec2(2:up_lim)...
			/abs(ref_H_rec2(round(1000/MIN_FREQ)+1)))),'b:',...
			x(2:up_lim+1),20*log10(abs(ref_H_TT(2:up_lim+1)/...
            abs(ref_H_TT(round(1000/MIN_FREQ)+1)))),'k:',...
			x(2:up_lim),20*log10(abs(H_mic(2:up_lim)))+ ...
			mic_level_correction,'g',...
			x(2:up_lim),20*log10(abs(H_rec1(2:up_lim)/...
			abs(H_rec1(round(1000/MIN_FREQ)+1)))),'r',...
			x(2:up_lim),20*log10(abs(H_rec2(2:up_lim)/...
			abs(H_rec2(round(1000/MIN_FREQ)+1)))),'b',...
			x(2:up_lim+1),20*log10(abs(H_TT(2:up_lim+1)/...
            abs(H_TT(round(1000/MIN_FREQ)+1)))),'k-'),
		axis([1 10000 -80 20]),grid on,
		text(15,17,strcat('green = mic, red = rec1, blue = rec2,',...
			' black=TT, ---- =current, ....=reference)')),
		zoom on,

		% plot absolute calibration (SNR & THD)
		h = findobj('Tag','snr_thd_fig');
		if (isempty(h)),
			snr_thd_fig,
		else,
			figure(h),
		end,
		semilogx(x(2:up_lim),ref_snr_rec1(2:up_lim),'rx', ...
			x(2:up_lim),ref_snr_rec2(2:up_lim),'bx',...
			x(2:up_lim),ref_snr_TT(2:up_lim),'kx',...
			x(2:up_lim),snr_rec1(2:up_lim),'ro', ...
			x(2:up_lim),snr_rec2(2:up_lim),'bo', ...
			x(2:up_lim),snr_TT(2:up_lim),'ko'),
		axis([1 10000 -50 140]), grid,
		text(12,115,['                red = rec1, ',...
			'  blue = rec2,  black = TT,  (o = current,  x= reference)']),
		text(90,105,['* mic  correction: ',...
			num2str(mic_level_correction,3),...
			'dB          (reference: ',...
			num2str(ref_mic_level_correction,3),'dB)']),
		text(90,95,['* max. level rec1: ',num2str(max_level_rec1,3),...
			'dB          (reference:  ',...
			num2str(ref_max_level_rec1,3),'dB)']),
		text(90,85,['* max. level rec2: ',num2str(max_level_rec2,3),...
			'dB          (reference:  ',...
			num2str(ref_max_level_rec2,3) 'dB)']),
		text(90,75,['* max. level TT:  ',num2str(max_level_TT,4),...
			'dB         (reference: ',...
			num2str(ref_max_level_TT,4),'dB)']),
		zoom on,


	case 'ok',
		try,
			load([OAE_PATH,'\DP_TT\DP_TT_calibration_data\00000__NONE__.mat']),
		catch,
			msgbox ('Do first calibrate!',...
                'ERROR calibration.m:','error'),
            return,
		end,
        
        save([OAE_PATH,...
                '\DP_TT\DP_TT_calibration_data\00000_last_calibration']...
            ,'H_rec1','H_rec2','H_TT','H_mic','snr_rec1','snr_rec2', ...
            'snr_TT','mic_level_correction','max_level_rec1','comment',...
            'max_level_rec2','max_level_TT','min_freq','sample_rate',...
			'BK_amplitude_94dBspl'),
        
        if (~exist([OAE_PATH,'\DP_TT\DP_TT_calibration_data\all']))
            eval(['!mkdir ', OAE_PATH,...
                    '\DP_TT\DP_TT_calibration_data\all'])
		end,
        
        save([OAE_PATH,'\DP_TT\DP_TT_calibration_data\all\cal',...
                num2str(now*10^9)],...
            'H_rec1','H_rec2','H_TT','H_mic','snr_rec1','snr_rec2',...
			'snr_TT','mic_level_correction','max_level_rec1','comment',...
			'max_level_rec2','max_level_TT','min_freq','sample_rate',...
			'BK_amplitude_94dBspl'),

		dp_tt_calibration('save'),
		
		
	case 'save', % save calibration as special name to keep as reference

		new_reference = get(findobj('Tag','new_name'),'String');
        if isempty(new_reference),
            return,
        end,
		
		try,
            load ([OAE_PATH,...
				'\DP_TT\DP_TT_calibration_data\00000__NONE__.mat']),
		catch,
			msgbox ('Do first calibrate!',...
                'ERROR calibration.m:','error'),
            return,
		end,
		
		str = [OAE_PATH,...
				'\DP_TT\DP_TT_calibration_data\', new_reference];
		if exist(str,'file'),
			if ~strcmp(questdlg(strvcat('File: ', str, ...
					' already exists. Save anyway?'),...
					'WARNING - dp_tt calibration:','No'),'Yes'),
				return,
			end,
		end
		try,
			save(str,'H_rec1','H_rec2','H_TT','H_mic','snr_rec1', ...
				'snr_rec2','snr_TT','mic_level_correction', ...
				'max_level_rec1','max_level_rec2','max_level_TT',...
				'min_freq','sample_rate','comment',...
			    'BK_amplitude_94dBspl'),
		catch,
			msgbox ('Reference name must be valid File name !',...
				'ERROR dp_tt calibration.m:','error'),
			return,
		end,
		dp_tt_calibration('fill'),
		

	case 'comment' % Callback of 'comment'-editbox of reference
        comment = get(findobj(gcbf,'Tag','comment'),'String');
        if ~iscell(comment),
            comment = comment;
        end,
        save([OAE_PATH,'\DP_TT\DP_TT_calibration_data\00000__NONE__.mat'],...
            'comment','-append'),
		
	case 'ref_comment' % Callback of 'comment'-editbox of reference
		h = findobj(gcbf,'Tag','references');
		str = get(h, 'String');
        if (isempty(str)), return, end,
        val = get(h, 'Value');
        comment = get(findobj(gcbf,'Tag','ref_comment'),'String');
        if ~iscell(comment),
            comment = comment;
        end,
        save([OAE_PATH,'\DP_TT\DP_TT_calibration_data\',str{val}],...
            'comment','-append'),
        
        
    case 'delete_reference',
        h = findobj('Tag','references');
        str = get(h, 'String');
        if (isempty(str)), return, end,
        val = get(h, 'Value');
        item = str{val};                
        if (strcmp(questdlg(strvcat('Are you sure you want delete : ',...
                item),'WARNING - DP_TT Calibration:','No'),'Yes')),
            delete ([OAE_PATH,'\DP_TT\DP_TT_calibration_data\',item])
            set(h,'Value',1);
            dp_tt_calibration('fill'),
		end 
		
        
    case 'close'
        h = findobj('Tag','phase_fig');
        if (~isempty(h)), close(h), end,
        h = findobj('Tag','magnitude_fig');
        if (~isempty(h)), close(h), end,
        h = findobj('Tag','snr_thd_fig');
        if (~isempty(h)), close(h), end,
        
        
    case 'help'
        msgbox (strvcat('ERROR DP_TT_calibration.m:',...
            'Not implemented yet. Sorry!')),
        
    otherwise
        error('Unkown action (case)!');
end %switch action     
    
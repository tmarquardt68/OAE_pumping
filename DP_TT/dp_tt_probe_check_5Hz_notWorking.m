%111111111111111111111111111111111111111111111111111111111111111111111112
% DP_TT_PROBE_CHECK(action) is part of the OAE toolbox (by Torsten Marquardt)
% and not of general use.
% Called by dp_tt_measurement.m

function dp_tt_probe_check(action)

global MIN_FREQ SAMPLE_RATE BITS_PER_SAMPLE OAE_PATH USE_CRIT
global CURRENT_EAR % will be filled by this function and is needed by 
% measure.m for the generation of korrektet primary and masker tones:
%  CURRENT_EAR(:,1) - rec1 transfer fct. (complex)
%  CURRENT_EAR(:,2) - rec2 transfer fct. (complex)
%  CURRENT_EAR(:,3) - TT   transfer fct. (complex)
%  CURRENT_EAR(:,4) - mic transfer fct. (complex)
%  CURRENT_EAR(1,5) - max. soundpressure level of rec1 [dB SPL]
%  CURRENT_EAR(2,5) - max. soundpressure level of rec2 [dB SPL]
%  CURRENT_EAR(3,5) - max. soundpressure level of TT [dB SPL]
%  CURRENT_EAR(:,6) - transfer fct.(complex) of current calibration


l = SAMPLE_RATE/MIN_FREQ; % length of a single sweep
n_avg = 5;

switch action
    case 'ini'
        dp_tt_probe_check_fig,
		drawnow,
        CURRENT_EAR(l,6) = 0;
        h = gcf;
        dp_tt_probe_check retry;
        uiwait(h),
        try,
			t=close(h);
		catch,
            CURRENT_EAR = [];
        end,
    
        
    case 'retry'
		USE_CRIT = 1;
		
        % load mat-file containing last calibration
        try,
			load([OAE_PATH ,...
                '\DP_TT\DP_TT_calibration_data\00000_last_calibration']);
		catch,
			CURRENT_EAR = [];
            msgbox(strvcat('Cannot find  data of current calibration: ',...
                'DP_TT\dp_tt_calibration_data\_last_calibration.mat!'...
                ,'New calibration required!'),'ERROR DP_TT_probe_check:',...
                'error');
            return,
        end,
        
        set(findobj('Tag','dp_tt_probe_check_fig'),...
            'Name','DP_TT Probe Check in progress (ca. 30 sec)')
        
        CURRENT_EAR(:,6) = H_rec1;
        CURRENT_EAR(:,4) = H_mic;
        CURRENT_EAR(1,5) = max_level_rec1;
        CURRENT_EAR(2,5) = max_level_rec2;
        CURRENT_EAR(3,5) = max_level_TT;
        
		% white noise generation
		rand('seed',0),
		theta = 2*pi*rand(1,l/2-1);
		[re im] = pol2cart([0 theta 0 ...
			-theta(l/2-1:-1:1)],... % "semi" pink noise
			ones(1,l)./[.1:0.1:l*.1]);
% 		[re im] = pol2cart([0 theta 0 ...
% 			-theta(l/2-1:-1:1)],...
% 			ones(1,l));
		H_wave = (re +im*i)';
		H_wave(1)=0;
% 		H_wave(round(10000/MIN_FREQ)+1:l-...
% 			round(10000/MIN_FREQ)+1)=0;
		
		H_wave_LF = H_wave;
		H_wave_LF(round(5/MIN_FREQ)+1:l-...
			round(5/MIN_FREQ)+1)=0;
		wave_LF = real(ifft(H_wave_LF));
		wave_LF = wave_LF./max(abs(wave_LF))/5;
        H_wave_LF = fft(wave_LF);

		wave = real(ifft(H_wave));
		wave = wave./max(abs(wave))/5;
		H_wave = fft(wave);
		
		% TF of receivers
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		% rec1
		[CURRENT_EAR(:,1), n_used_rec1] = data_acquisition(wave, 1, 1, ...
			ceil(n_avg*MIN_FREQ), USE_CRIT);
		CURRENT_EAR(:,1)= fft(CURRENT_EAR(:,1))./H_wave./H_mic;

        % rec2
        [CURRENT_EAR(:,2), n_used_rec2]  = data_acquisition(...
            [zeros(length(wave),1), wave], 2, 1,ceil(n_avg*MIN_FREQ), USE_CRIT);
		CURRENT_EAR(:,2)= fft(CURRENT_EAR(:,2))./H_wave./H_mic;
      
        % TT
        [CURRENT_EAR(:,3), n_used_TT] = data_acquisition(...
            [zeros(length(wave),2), wave_LF], ...
			3, 1, ceil(n_avg*MIN_FREQ), USE_CRIT);
		CURRENT_EAR(:,3)= fft(CURRENT_EAR(:,3))./H_wave_LF./H_mic;
				
		% plot magnitude relative to last calibration
		h_axes(1) = findobj('Tag','axes_magnitude');
		axes(h_axes(1)),
		zoom on,
		x = linspace(0,SAMPLE_RATE/2,l/2);
		semilogx(...
			x(5/MIN_FREQ+1:l/2),...
			20*log10(abs(CURRENT_EAR(5/MIN_FREQ+1:l...
			/2,1)./H_rec1(5/MIN_FREQ+1:l/2))),'b',...
			x(5/MIN_FREQ+1:l/2),...
            20*log10(abs(CURRENT_EAR(5/MIN_FREQ+1:l...
            /2,2)./H_rec2(5/MIN_FREQ+1:l/2))),'r',...
			x(5/MIN_FREQ+1:5/MIN_FREQ),...
            20*log10(abs(CURRENT_EAR(5/MIN_FREQ+1:5/MIN_FREQ,3)...
            ./H_TT(5/MIN_FREQ+1:5/MIN_FREQ))),'k',...
			x,zeros(1,length(x)),'g:'), 
        axis([10 24000 -20 20]),grid on
        text(20, 9, strcat('MAGNITUDE:    ',...
            'blue=rec1, red=rec2, black= TT, green = reference (all in dB)'),...
			'FontSize', 10),
        set(h_axes(1),'Tag','axes_magnitude'),
        
%         % plot phase
%         h_axes(2) = findobj('Tag','axes_phase');
%         axes(h_axes(2)),
%         semilogx(x(20/MIN_FREQ+1:10000/MIN_FREQ),unwrap(angle...
%             (CURRENT_EAR(20/MIN_FREQ+1:10000/MIN_FREQ,1))-...
%             angle(H_rec1(20/MIN_FREQ+1:10000/MIN_FREQ))),'b',...
%             x(20/MIN_FREQ+1:10000/MIN_FREQ),unwrap(angle...
%             (CURRENT_EAR(20/MIN_FREQ+1:10000/MIN_FREQ,2))-...
%             angle(H_rec2(20/MIN_FREQ+1:10000/MIN_FREQ))),'r'),
%         text(10,0,'PHASE as radiant'),
%         set(h_axes(2),'Tag','axes_phase'),

        % plot absolute magnitude 
        h_axes(2) = findobj('Tag','axes_phase');
        axes(h_axes(2)),
		semilogx(...
			x(5/MIN_FREQ+1:l/2),...
			20*log10(abs(CURRENT_EAR(5/MIN_FREQ+1:l...
			/2,1)./CURRENT_EAR(1000/MIN_FREQ+1,1)))+CURRENT_EAR(1,5),'b',...
			x(5/MIN_FREQ+1:l/2),...
            20*log10(abs(CURRENT_EAR(5/MIN_FREQ+1:l...
            /2,2)./CURRENT_EAR(1000/MIN_FREQ+1,2)))+CURRENT_EAR(2,5),'r',... 
			x(5/MIN_FREQ+1:l/2),...
            20*log10(abs(CURRENT_EAR(5/MIN_FREQ+1:l...
            /2,3)./CURRENT_EAR(5/MIN_FREQ+1,3)))+CURRENT_EAR(3,5),'k'), 
        grid on, zoom on,
        set(h_axes(2),'Tag','axes_phase'),
        linkaxes(h_axes,'x')        
        axis([10 24000 40 125]),grid on


        % normalize to _last_calibration at 1kHz / 20 Hz
        CURRENT_EAR(:,1) = CURRENT_EAR(:,1)./...
            abs(H_rec1(round(1000/MIN_FREQ)+1));
        CURRENT_EAR(:,2)=CURRENT_EAR(:,2)./...
            abs(H_rec2(round(1000/MIN_FREQ)+1));
        CURRENT_EAR(:,3)=CURRENT_EAR(:,3)./...
            abs(H_TT(round(5/MIN_FREQ)+1));
        
        set(findobj('Tag','dp_tt_probe_check_fig'),...
            'Name','DP_TT Probe Check finished. OK or RETRY?')
		
        
    case 'ok'
        save([OAE_PATH, '\DP_TT\_last_probe_check'], 'CURRENT_EAR'),
        set(gcf,'Tag','ok is pressed');
        uiresume,
        
    case 'load_last'
        % load mat-file containing last calibration
%         h = warndlg(strvcat('Are you sure the probe setting has not changed  ',...
%             'sincethe last probe fit check? The probe check data will´',...
%             'be used to normalize the mesurement results!',...
%         	'It is recommended to do always a probe fit check!'),...
%             'WARNING: DP_TT_probe_check:');
%         uiwait(h),
% 		drawnow,
        try,
			load([OAE_PATH ,'\DP_TT\_last_probe_check']),
		catch,
            h = errordlg(strvcat('Cannot find  data of last probe fit check: ',...
                'DP_TT\_last_probe_check.mat!'...
                ,'Do probe_check now!'),'ERROR DP_TT_probe_check:');
            uiwait(h),
                CURRENT_EAR = []; % empty means probe_check failed
            return,
        end,
    otherwise
        error('Unkown action (case)!');
end     % switch action        

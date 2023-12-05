% DP_TT_viewer(action, patient, measurement):
%
%     patient:        patient name
%     measurement:    file name of measurement
%        Both will be save in the figure property 'User' using
%        a cell array named 'series_spec' (series_spec{1} = patient,
%        series_spec{2} = measurement)
% DP_TT_VIEWER is part of the OAE toolbox (by Torsten Marquardt)
% and  not of general use.
% MODIFICATIONS: [98-08-04H=Hensel lines 89,92
% (noise_mask was different in 'plot_hold' vs 'plot_replace')];
% 26.09.98 T:uses right and left neighbour of each modulation line
%            to calculate time course of noise ([right(t) + left(t)]/2)
% FURTHER MODIFICATIONS:98-10-05H LINE381,98-11-12LINES18,189-192,381-385
function dp_tt_viewer(action, patient, measurement)

global OAE_PATH
henselnoise=0; %MODIFIED98-11-12 display of Noise Floor in LINES189-192,381-385

switch action
    case 'ini'
        dp_tt_viewer_fig,
        series_spec{1} = patient;
        series_spec{2} = measurement;
        set(gcf,'User',series_spec);
        set(gcf,'Name',['DP_TT_Viewer: ', patient,' , ' measurement(6:20)]),
        dp_tt_viewer('fill',patient,measurement),


    case 'fill'
        h_fig = findobj('Name',['DP_TT_Viewer: ', patient,' , ',measurement(6:20)]);
        series_spec = get(h_fig,'User');
        figure(h_fig),
        load([OAE_PATH,'\Subjects\',patient,'\',measurement],'F1',...
            'L1','l1','F2','L2','l2','F_TT','L_TT','l_TT','l_dp','phase_dp',...
            'comment','min_freq'),
        no_phase_dp = ~exist('phase_dp','var');
        n = length(l1);
        if (~exist('l_dp','var'))
            l_dp=-99*ones(1,n);
        elseif size(l_dp,2) == 1
            l_dp = [l_dp' zeros(n,4)];
        end,
        set(findobj(h_fig,'Tag','comment'),'String',comment),
        if (isempty(l1)),
            return,
        end,

        % plot  series
        string = ''; n_dp=0; n_sf=0;series_DP = [];series_CM = [];series_SF = [];
        for n = 1:length(l1(~isnan(l1)))


            if strcmp(get(gcbo,'String'),'View') %&& no_phase_dp
                load([OAE_PATH,'\Subjects\',series_spec{1},'\', series_spec{2}], ...
                    ['avg', num2str(n)], 'h_mic'),
                eval(['avg = avg',num2str(n) ';']);
                H_avg_acoust = fft(avg(:,1))./h_mic;
                if size(avg,2)>1
                    H_avg_CM = fft(avg(:,2));
                else
                    H_avg_CM = zeros(size(avg,1),1);
                end

                % get the actual levels
                if L1(n) < 0 % SFOAE
                    f(1) = F2(n)+1;
                    f(2) = 0;
                    f(3) = F2(n)-2*F_TT(n)+1;
                    f(4) = F2(n)-F_TT(n)+1;
                    f(5) = F2(n)+F_TT(n)+1;
                    f(6) = F2(n)+2*F_TT(n)+1;
                    f(7) = 0;
                    f(8) = 0;
                    f(9) = 0;
                    f(10) = 0;
                else % DPOAE
                    f(1) = 2*F1(n)-F2(n)+1;
                    f(2) = F2(n)-F1(n)+1;
                    f(3) = 2*F1(n)-F2(n)-2*F_TT(n)+1;
                    f(4) = 2*F1(n)-F2(n)-F_TT(n)+1;
                    f(5) = 2*F1(n)-F2(n)+F_TT(n)+1;
                    f(6) = 2*F1(n)-F2(n)+2*F_TT(n)+1;
                    f(7) = F2(n)-F1(n)-2*F_TT(n)+1;
                    f(8) = F2(n)-F1(n)-F_TT(n)+1;
                    f(9) = F2(n)-F1(n)+F_TT(n)+1;
                    f(10) = F2(n)-F1(n)+2*F_TT(n)+1;
                end
                for q=length(f):-1:1
                    if f(q)>0 && f(q) < length (H_avg_acoust)
                        l_dp(n,q) = round(20*log10(abs(H_avg_acoust(f(q))))*10)/10;
                        phase_dp(n,q) = unwrap(angle(H_avg_acoust(f(q))))/2/pi;
                        l_dp(n,q+10) = round(20*log10(abs(H_avg_CM(f(q))))*10)/10;
                        phase_dp(n,q+10) = unwrap(angle(H_avg_CM(f(q))))/2/pi;
                    else
                        l_dp(n,q) = nan;
                        l_dp(n,q+10) = nan;
                        phase_dp(n,q) = nan;
                        phase_dp(n,q+10) = nan;
                    end
                end
                l_TT(n) = round(20*log10(abs(H_avg_acoust(F_TT(n)+1)))*10)/10;
                l1(n) = round(20*log10(abs(H_avg_acoust(F1(n)+1))));
                l2(n) = round(20*log10(abs(H_avg_acoust(F2(n)+1))));
            end


            
            str_TT = sprintf(' %5.1f',l_TT(n));
            if (abs(l_TT(n)-L_TT(n)) > 0),
                str_TT = [str_TT,sprintf('(%3d)',L_TT(n))];
            else,
                str_TT = [str_TT,'     '];
            end;

            str_L1 = sprintf(' %2d',l1(n));
            if (abs(l1(n)-L1(n)) > 0),
                str_L1 = [str_L1,sprintf('(%3d)',L1(n))];
            else,
                str_L1 = [str_L1,'     '];
            end;

            str_L2 = sprintf(' %2d',l2(n));
            if (abs(l2(n)-L2(n)) > 0),
                str_L2 = [str_L2,sprintf('(%2d)',L2(n))];
            else,
                str_L2 = [str_L2,'    '];
            end;

            if L1(n)>= 0 % DPAOE
                n_dp = n_dp+1;
                series_DP(n_dp,:) = [n, floor(F_TT(n)*min_freq),l_TT(n),L_TT(n),floor(F1(n)*min_freq),...
                    l1(n),L1(n),floor(F2(n)*min_freq),l2(n),L2(n),l_dp(n,:),phase_dp(n,:)];
            else % SFOAE
                n_sf = n_sf+1;
                series_SF(n_sf,:) = [n, floor(F_TT(n)*min_freq),l_TT(n),L_TT(n),floor(F1(n)*min_freq),...
                    l1(n),L1(n),floor(F2(n)*min_freq),l2(n),L2(n),l_dp(n,:),phase_dp(n,:)];
            end

            string{n}=sprintf(....
                '%2d%3d%1s%5d%1s%5d%1s%5.1f%5.1f%5.1f%5.1f%5.1f%5.1f',...
                n, floor(F_TT(n)*min_freq),str_TT,floor(F1(n)*min_freq),...
                str_L1,floor(F2(n)*min_freq),str_L2,l_dp(n,1),l_dp(n,2),...
                l_dp(n,3),l_dp(n,4),l_dp(n,5),l_dp(n,6));
        end
        save([OAE_PATH,'\Subjects\',series_spec{1},'\', series_spec{2}], ...
            'l_dp','phase_dp', '-append'),

        %% Plot series
        position = [1920   20   560   1100];
        h_series = findobj('Position',position);
        if (isempty(h_series))
            h_series = figure;
        end
        figure(h_series), clf,
        set(h_series,'Position',position,...
            'Name',['Series: ',series_spec{1}(1:5),...
            ' , ',series_spec{2}(6:20)])
        h_ax = []; a=1;
        if ~isempty(series_SF)
            if size(l_dp,2)>10 % modCM are recorded
                % Plot modCM level
                h_ax(a)=subplot('position',[0.0446    0.1443    0.9429    0.0984]); a=a+1;
                plot(series_SF(:,23),'bx-'), hold on
                plot(series_SF(:,24),'rx-')
                plot(series_SF(:,21),'kx-','LineWidth',2),
                plot(series_SF(:,25),'rx-')
                plot(series_SF(:,26),'bx-')
                grid on,title modCM,zoom on
                % Plot modCM phase
                h_ax(a)=subplot('position',[0.0446    0.0318    0.9429    0.0984]); a=a+1;
                plot(series_SF(:,23+20),'bx-'), hold on
                plot(series_SF(:,24+20),'rx-')
                plot(series_SF(:,21+20),'kx-','LineWidth',2),
                plot(series_SF(:,25+20),'rx-')
                plot(series_SF(:,26+20),'bx-')
                grid on,zoom on
                xlabel('no of 5s-presentions'),zoom on
                % Plot modSF level
                h_ax(a)=subplot('position',[0.0446    0.3902    0.9429    0.0984]); a=a+1;
                hold off
                plot(series_SF(:,13),'bx-')
                hold on
                plot(series_SF(:,14),'rx-')
                plot(series_SF(:,15),'rx-')
                plot(series_SF(:,16),'bx-'),
                grid on,title modSF,zoom on
                % Plot modSF phase
                h_ax(a)=subplot('position',[0.0446    0.2741    0.9429    0.0984]); a=a+1;
                hold off
                plot(series_SF(:,13+20),'bx-')
                hold on
                plot(series_SF(:,14+20),'rx-')
                plot(series_SF(:,15+20),'rx-')
                plot(series_SF(:,16+20),'bx-'),
                grid on,zoom on
            else % plot only mod_SF (CM not recored)
                % Plot modSF level
                h_ax(a)=subplot('position',[0.0446    0.3902    0.9429    0.0984]); a=a+1;
                hold off
                plot(series_SF(:,13),'bx-')
                hold on
                plot(series_SF(:,14),'rx-')
                plot(series_SF(:,15),'rx-')
                plot(series_SF(:,16),'bx-'),
                grid on,title modSF,zoom on
                % Plot modSF phase
                h_ax(a)=subplot('position',[0.0446    0.2741    0.9429    0.0984]); a=a+1;
                hold off
                plot(series_SF(:,13+10),'bx-')
                hold on
                plot(series_SF(:,14+10),'rx-')
                plot(series_SF(:,15+10),'rx-')
                plot(series_SF(:,16+10),'bx-'),
                grid on,zoom on
                subplot(6,1,5,'replace')
                subplot(6,1,6,'replace')
                xlabel('no of 5s-presentions'),zoom on
            end
        else
            subplot('position',[0.0446    0.3902    0.9429    0.0984],'replace')
            subplot('position',[0.0446    0.2741    0.9429    0.0984],'replace')
            subplot('position',[0.0446    0.1443    0.9429    0.0984],'replace')
            subplot('position',[0.0446    0.0318    0.9429    0.0984],'replace')
            xlabel('no of 5s-presentions'),zoom on
        end
        if ~isempty(series_DP)
            if size(l_dp,2)>10  % modCM, but no modSF recorded
                if isempty(series_SF)
                    % Plot modCM level
                    h_ax(a)=subplot('position',[0.0446    0.1443    0.9429    0.0984]); a=a+1;
                    plot(series_DP(:,23),'bx-'), hold on
                    plot(series_DP(:,24),'rx-')
                    plot(series_DP(:,21),'kx-','LineWidth',2),
                    plot(series_DP(:,25),'rx-')
                    plot(series_DP(:,26),'bx-')
                    grid on,title modCM,zoom on
                    % Plot modCM phase
                    h_ax(a)=subplot('position',[0.0446    0.0318    0.9429    0.0984]); a=a+1;
                    plot(series_DP(:,23+20),'bx-'), hold on
                    plot(series_DP(:,24+20),'rx-')
                    plot(series_DP(:,21+20),'kx-','LineWidth',2),
                    plot(series_DP(:,25+20),'rx-')
                    plot(series_DP(:,26+20),'bx-')
                    grid on,zoom on
                    xlabel('no of 5s-presentions'),zoom on
                end
                % Plot modDP level
                h_ax(a)=subplot('position',[0.0446    0.7527    0.9429    0.2337]); a=a+1;
                hold off
                plot(series_DP(:,11),'kx-', 'LineWidth', 2),
                hold on
                plot(series_DP(:,13),'bx-')
                plot(series_DP(:,14),'rx-')
                plot(series_DP(:,15),'rx-')
                plot(series_DP(:,16),'bx-')
                plot(series_DP(:,12),'ko:', 'LineWidth', 2),
                plot(series_DP(:,17),'bo:')
                plot(series_DP(:,18),'ro:')
                plot(series_DP(:,19),'ro:')
                plot(series_DP(:,20),'bo:')
                grid on,title modDP,zoom on
                if size(series_DP,1)>2
                    ylim([min(series_DP(2:end-1,13:20),[],'all') max(series_DP(:,11:20),[],'all')])
                end
                    % Plot modDP phase
                h_ax(a)=subplot('position',[0.0446    0.5164    0.9429    0.2191]); a=a+1;
                hold off
                plot(series_DP(:,11+20),'kx-', 'LineWidth', 2),
                hold on
                plot(series_DP(:,13+20),'bx-')
                plot(series_DP(:,14+20),'rx-')
                plot(series_DP(:,15+20),'rx-')
                plot(series_DP(:,16+20),'bx-')
                plot(series_DP(:,12+20),'ko:', 'LineWidth', 2),
                plot(series_DP(:,17+20),'bo:')
                plot(series_DP(:,18+20),'ro:')
                plot(series_DP(:,19+20),'ro:')
                plot(series_DP(:,20+20),'bo:')
                grid on,zoom on
            else % plot only modDP (CM not recored)
                % Plot modDP level
                h_ax(a)=subplot('position',[0.0446    0.7527    0.9429    0.2337]); a=a+1;
                hold off
                plot(series_DP(:,11),'kx-', 'LineWidth', 2),
                hold on
                plot(series_DP(:,13),'bx-')
                plot(series_DP(:,14),'rx-')
                plot(series_DP(:,15),'rx-')
                plot(series_DP(:,16),'bx-')
                plot(series_DP(:,12),'ko:', 'LineWidth', 2),
                plot(series_DP(:,17),'bo:')
                plot(series_DP(:,18),'ro:')
                plot(series_DP(:,19),'ro:')
                plot(series_DP(:,20),'bo:')
                grid on,title modDP,zoom on
                ylim([min(series_DP(2:end-1,13:20),[],'all') max(series_DP(:,11:20),[],'all')])
                % Plot modDP phase
                h_ax(a)=subplot('position',[0.0446    0.5164    0.9429    0.2191]); a=a+1;
                hold off
                plot(series_DP(:,11+10),'kx-', 'LineWidth', 2),
                hold on
                plot(series_DP(:,13+10),'bx-')
                plot(series_DP(:,14+10),'rx-')
                plot(series_DP(:,15+10),'rx-')
                plot(series_DP(:,16+10),'bx-')
                plot(series_DP(:,12+10),'ko:', 'LineWidth', 2),
                plot(series_DP(:,17+10),'bo:')
                plot(series_DP(:,18+10),'ro:')
                plot(series_DP(:,19+10),'ro:')
                plot(series_DP(:,20+10),'bo:')
                grid on,zoom on
            end
        else
            subplot('position',[0.0446    0.7527    0.9429    0.2337],'replace')
            subplot('position',[0.0446    0.5164    0.9429    0.2191],'replace')
        end
        h = findobj(h_fig,'Tag','atomics');
        set(h,'String',string),
        set(h,'Value',n);
        linkaxes(h_ax, 'x')

    case 'plot_replace'
        if exist('patient','var')
            h_fig = findobj('Name',['DP_TT_Viewer: ', patient,' , ',measurement(6:20)]);
        else
            h_fig = gcbf;
        end
        series_spec = get(h_fig,'User');
        scrsz = get(0,'ScreenSize');
        n = get(findobj(h_fig,'Tag','atomics'),'Value');
        % get order of modulation
        order = str2num(get(findobj(h_fig,'Tag','order'),'String'));
        if isempty(order)
            msgbox('define the order of modulation',...
                'ERROR: DP_TT Viewer','error')
            return,
        end

        load([OAE_PATH,'\Subjects\',series_spec{1},'\',series_spec{2}],...
            'min_freq', 'F1','F2','F_TT','L_TT','L1','l_dp'),

        load([OAE_PATH,'\Subjects\',series_spec{1},'\', series_spec{2}], ...
            ['avg', num2str(n)], 'h_mic'),
        eval(['H_atomic = fft(avg',num2str(n), '(:,1))./h_mic;']);
        try
        eval(['H_eCoch = fft(avg',num2str(n), '(:,2));']);
        catch
        end
        % Use for old measurments of eg. Manuela Berger und Jana Wolf
        % load([OAE_PATH,'\Subjects\',series_spec{1},'\', series_spec{2}],'gains'),
        % eval(['H_atomic = (fft(-avg',num2str(n),'./10^(gains(n)/20))./h_mic)'';']);
        % AND REPLACE line 392
        %    course = abs(ifft(H_atomic .* mod_mask));
        %    course = abs(ifft((H_atomic .* mod_mask)'));
        % i.e. IFFT gives reversed time course when applying to transposed matrix!
        % END: Use for old measurments of eg. Manuela Berger und Jana Wolf

        l = length(H_atomic);
        noise_mask1(l,1) = 0;
        noise_mask2(l,1) = 0;

        if L1(n)>= 0 % plot DPOAE spectrum
            mod_mask(l,1) = 0;
            f(1) = 2*F1(n)-F2(n)-2*F_TT(n)+1;
            f(2) = 2*F1(n)-F2(n)-F_TT(n)+1;
            f(3) = 2*F1(n)-F2(n)+F_TT(n)+1;
            f(4) = 2*F1(n)-F2(n)+2*F_TT(n)+1;
            f(5) = F2(n)-F1(n)-2*F_TT(n)+1;
            f(6) = F2(n)-F1(n)-F_TT(n)+1;
            f(7) = F2(n)-F1(n)+F_TT(n)+1;
            f(8) = F2(n)-F1(n)+2*F_TT(n)+1;
            f(9) = F1(n)-2*F_TT(n)+1;
            f(10) = F1(n)-F_TT(n)+1;
            f(11) = F1(n)+F_TT(n)+1;
            f(12) = F1(n)+2*F_TT(n)+1;
            f(13) = F2(n)-2*F_TT(n)+1;
            f(14) = F2(n)-F_TT(n)+1;
            f(15) = F2(n)+F_TT(n)+1;
            f(16) = F2(n)+2*F_TT(n)+1;
            for q=1:length(f)
                if f(q)>0 & f(q)< l
                    mod_mask(f(q)) = 1;
                end
            end

            position = [0,round(scrsz(4)*3/4)-18,round(scrsz(3)/3.6),...
                round(scrsz(4)/4-40)];
            h_spectrum = findobj('Position',position);
            if (isempty(h_spectrum)),
                h_spectrum = viewer_spctr_fig;
            end,
            figure(h_spectrum), clf,
            x = [0:min_freq:l*min_freq];
            hold on,

            stem_60(x(1:end-1),20*log10(abs(H_atomic(1:length(x)-1))),'kx'),
            stem_60(x(1:F_TT(n):length(x)-1),20*log10(abs(H_atomic(1:F_TT(n):length(x)-1))),'y.'),

            stem_60(x([F_TT(n)+1,F1(n)+1,F2(n)+1,2*F1(n)-F2(n)+1,F2(n)-F1(n)+1]),20*log10(abs(H_atomic...
                ([F_TT(n)+1,F1(n)+1,F2(n)+1,2*F1(n)-F2(n)+1,F2(n)-F1(n)+1]))),'filled','bo'),

            if 2*F1(n)-F2(n)+1 > 0 && 2*F1(n)-F2(n)+1 < length(x)-1
                 stem_60(x(2*F1(n)-F2(n)+1),20*log10(abs(H_atomic(2*F1(n)-F2(n)+1))),'filled','bo'),
            end
            if F2(n)-F1(n)+1 > 0 && F2(n)-F1(n)+1 < length(x)-1
                 stem_60(x(F2(n)-F1(n)+1),20*log10(abs(H_atomic(F2(n)-F1(n)+1))),'filled','bo'),
            end

            stem_60(x(mod_mask>0),20*log10(abs(H_atomic(mod_mask>0))),'filled','ro'),

            harmonics_mask = []; q=2;
            while q*F1(n) < length(x)-1
                harmonics_mask = [harmonics_mask q*F1(n)];
                q=q+1;
            end
            stem_60(x(harmonics_mask+1),20*log10(abs(H_atomic...
                (harmonics_mask+1))),'b.'),
            harmonics_mask = []; q=2;
            while q*F2(n) < length(x)-1
                harmonics_mask = [harmonics_mask q*F2(n)];
                q=q+1;
            end
            stem_60(x(harmonics_mask+1),20*log10(abs(H_atomic...
                (harmonics_mask+1))),'b.'),
            stem_60(x(1:F_TT(n):length(x)-1),20*log10(abs(H_atomic(1:F_TT(n):length(x)-1))),'y.'),

            axis([0 24000 -40 115]),grid on, zoom on,
            set(h_spectrum,'Position',position,...
                'Name',['DP_TT_spctr: Nr.',num2str(n),', ',series_spec{1},...
                ' , ',series_spec{2}(6:20)]),

            % plot CM spectrum
            position = [round(scrsz(3)/3.6),round(scrsz(4)*3/4)-18,...
                round(scrsz(3)/3.6),round(scrsz(4)/4-40)];
            h_spectrum = findobj('Position',position);
            if (isempty(h_spectrum)),
                h_spectrum = viewer_spctr_fig;
            end,
            figure(h_spectrum), clf,
            x = [0:min_freq:l*min_freq];
            hold on,
            offset = 0;
            stem_60(x(1:end-1),offset + 20*log10(abs(H_eCoch(1:length(x)-1,1))),'kx'),
            stem_60(x([F_TT(n)+1,F1(n)+1,F2(n)+1,2*F1(n)-F2(n)+1,F2(n)-F1(n)+1]),...
                offset + 20*log10(abs(H_eCoch...
                ([F_TT(n)+1,F1(n)+1,F2(n)+1,2*F1(n)-F2(n)+1,F2(n)-F1(n)+1]))),'filled','bo'),
            stem_60(x(mod_mask>0),offset + 20*log10(abs(H_eCoch(mod_mask>0))),'filled','ro'),
            harmonics_mask = []; q=2;
            while q*F2(n) < length(x)-1
                harmonics_mask = [harmonics_mask q*F2(n)];
                q=q+1;
            end
            stem_60(x(harmonics_mask+1),20*log10(abs(H_eCoch(harmonics_mask+1))),'b.'),
            stem_60(x(1:F_TT(n):length(x)-1),20*log10(abs(H_eCoch(1:F_TT(n):length(x)-1))),'y.'),
          
            axis([0 24000 -40 100]),grid on, zoom on,
            set(h_spectrum,'Position',position,...
                'Name',['CM_spctr: Nr.',num2str(n),', ',series_spec{1},...
                ' , ',series_spec{2}(6:20)]),

        else % plot SFOAE spectrum 
            mod_mask(l,1) = 0;
            f(1) = F2(n)-2*F_TT(n)+1;
            f(2) = F2(n)-F_TT(n)+1;
            f(3) = F2(n)+F_TT(n)+1;
            f(4) = F2(n)+2*F_TT(n)+1;
            for q=1:4
                if f(q)>0 && f(q)< l
                    mod_mask(f(q)) = 1;
                end
            end

            position = [0,round(scrsz(4)*3/4)-18,round(scrsz(3)/3.6),...
                round(scrsz(4)/4-40)];
            h_spectrum = findobj('Position',position);
            if (isempty(h_spectrum)),
                h_spectrum = viewer_spctr_fig;
            end,
            figure(h_spectrum), clf,
            x = [0:min_freq:l*min_freq];
            hold on,

            stem_60(x(1:end-1),20*log10(abs(H_atomic(1:length(x)-1,1))),'kx'),
            % stem_60(x(1:F_TT(n):length(x)),20*log10(abs(H_atomic(1:F_TT(n):length(x)))),'.'),

            stem_60(x([F_TT(n)+1,F2(n)+1]),20*log10(abs(H_atomic...
                ([F_TT(n)+1,F2(n)+1]))),'filled','bo'),

            stem_60(x(mod_mask>0),20*log10(abs(H_atomic(mod_mask>0))),'filled','ro'),

            harmonics_mask = []; q=2;
            while q*F2(n) < length(x)-1
                harmonics_mask = [harmonics_mask q*F2(n)];
                q=q+1;
            end
            stem_60(x(harmonics_mask+1),20*log10(abs(H_atomic(harmonics_mask+1))),'b.'),
            axis([0 24000 -40 115]),grid on, zoom on,
            set(h_spectrum,'Position',position,...
                'Name',['DP_TT_spctr: Nr.',num2str(n),', ',series_spec{1},...
                ' , ',series_spec{2}(6:20)]),

            % plot CM spectrum
            position = [round(scrsz(3)/3.6),round(scrsz(4)*3/4)-18,...
                round(scrsz(3)/3.6),round(scrsz(4)/4-40)];
            h_spectrum = findobj('Position',position);
            if (isempty(h_spectrum)),
                h_spectrum = viewer_spctr_fig;
            end,
            figure(h_spectrum), clf,
            x = [0:min_freq:l*min_freq];
            hold on,
            offset = 0;
            stem_60(x(1:end-1),offset + 20*log10(abs(H_eCoch(1:length(x)-1,1))),'kx'),
            stem_60(x([F_TT(n)+1,F2(n)+1]),offset + 20*log10(abs(H_eCoch...
                ([F_TT(n)+1,F2(n)+1]))),'filled','bo'),
            stem_60(x(mod_mask>0),offset + 20*log10(abs(H_eCoch(mod_mask>0))),'filled','ro'),
            harmonics_mask = []; q=2;
            while q*F2(n) < length(x)-1
                harmonics_mask = [harmonics_mask q*F2(n)];
                q=q+1;
            end
            stem_60(x(harmonics_mask+1),20*log10(abs(H_eCoch(harmonics_mask+1))),'b.'),
            stem_60(x(1:F_TT(n):length(x)-1),20*log10(abs(H_eCoch(1:F_TT(n):length(x)-1))),'y.'),
          
            axis([0 24000 -40 100]),grid on, zoom on,
            set(h_spectrum,'Position',position,...
                'Name',['CM_spctr: Nr.',num2str(n),', ',series_spec{1},...
                ' , ',series_spec{2}(6:20)]),
        end


%         % plot magnitude
%         position = [round(scrsz(3)/3.6*2),round(scrsz(4)*3/4)-18,...
%         round(scrsz(3)/3.6),round(scrsz(4)/4-40)];
%         h_magnitude = findobj('Position',position);
%         if (isempty(h_magnitude)),
%             h_magnitude = viewer_magn_fig;
%         end,
%         figure(h_magnitude), hold off,
%         x = linspace(0,2/F_TT(n)/min_freq*1000,round(2*l/F_TT(n)));
%         if L1(n)>= 0 % % plot time course of DPOAE magintude
%             hold off
%             mod_mask = zeros(l,1);
%             mod_mask([...
%                 2*F1(n)-F2(n)+1, ...
%                 2*F1(n)-F2(n)+1-2*F_TT(n), ...
%                 2*F1(n)-F2(n)+1-F_TT(n), ...
%                 2*F1(n)-F2(n)+1+F_TT(n), ...
%                 2*F1(n)-F2(n)+1+2*F_TT(n), ...
%                 ]) = 1;
%             course = abs(ifft((H_atomic .* mod_mask)));
%             plot(x,20*log10(l * course(1:round(2*l/F_TT(n))))),
%             hold on;
%             try
%                 mod_mask = zeros(l,1);
%                 mod_mask([...
%                     F2(n)-F1(n)+1, ...
%                     F2(n)-F1(n)+1-2*F_TT(n), ...
%                     F2(n)-F1(n)+1-F_TT(n), ...
%                     F2(n)-F1(n)+1+F_TT(n), ...
%                     F2(n)-F1(n)+1+2*F_TT(n), ...
%                     ]) = 1;
%                 course = abs(ifft((H_atomic .* mod_mask)));
%                 plot(x,20*log10(l * course(1:round(2*l/F_TT(n)))),'b--'),
%                 legend({'2f2-f1' 'f2-f1'})
%                 xlim([0 2/F_TT(n)/min_freq*1000]), grid on, zoom on,
%                 set(h_magnitude,'Position',position,...
%                     'Name',['DP_TT_mag: Nr.',num2str(n),', ',series_spec{1},' , ',...
%                     series_spec{2}(6:20)]),
%             catch
%             end
%         else % time course of CM
%             hold off
%             try 
%             mod_mask = zeros(l,1);
%             mod_mask([...
%                 F2(n)+1, ...
%                 F2(n)+1-2*F_TT(n), ...
%                 F2(n)+1-F_TT(n), ...
%                 F2(n)+1+F_TT(n), ...
%                 F2(n)+1+2*F_TT(n), ...
%                 ]) = 1;
%             course = abs(ifft((H_atomic .* mod_mask)));
%             plot(x,20*log10(l * course(1:round(2*l/F_TT(n))))),
%             xlim([0 2/F_TT(n)/min_freq*1000]), grid on, zoom on,
%             set(h_magnitude,'Position',position,...
%                 'Name',['CM_mag: Nr.',num2str(n),', ',series_spec{1},' , ',...
%                 series_spec{2}(6:20)]),
%             catch
%             end
%       end
%         if(L_TT(n)>50)
%             H_f_TT(l) = 0;
%             H_f_TT(F_TT(n)+1) = H_atomic(F_TT(n)+1);
%             course_TT = real(ifft(H_f_TT));
%             course_TT = course_TT/max(course_TT);
%             yLimit = get(gca,'YLim');
%             hold on
%             plot(x, mean(yLimit) + diff(yLimit)/2 * course_TT(1:round(2*l/F_TT(n))),'k:'),
%         end


        %   % plot phase
        %     position = [round(scrsz(3)/3.6*2),round(scrsz(4)*3/4)-18,...
        %             round(scrsz(3)/3.6),round(scrsz(4)/4-40)];
        %     h_phase = findobj('Position',position);
        %     if (isempty(h_phase)),
        %         h_phase = viewer_phase_fig,
        % 	end,
        %     figure(h_phase), hold off,
        %     f_DP_mask = zeros(l,1);
        %     f_DP_mask(f_DP + 1) = 1;
        %     course = unwrap(angle(ifft(H_atomic.*mod_mask)))...
        %         - unwrap(angle(ifft(H_atomic.*f_DP_mask)));
        %     x = linspace(0,1/F_TT(n)/min_freq*1000,round(l/F_TT(n)));
        %     % plot time course of DP phase
        %     plot(x,(course(1:round(l/F_TT(n))))),
        %     hold on,zoom on,
        %     % plot time course of the masker
        %     if(L_TT(n)>50),
        %             Y = get(gca,'YTick');
        %             tmp = Y(ceil(length(Y)/2));
        % 			plot(x,(tmp - Y(1))*course_TT(1:round(l/F_TT(n)))+tmp,'k:'),
        % 	end,
        % 	grid on;
        % 	set(h_phase,'Position',position, ...
        % 		'Name',['DP_TT_phase: Nr.',num2str(n),', Ord. ',...
        % 		num2str(order),', ',series_spec{1}(1:3),' , ',...
        % 		datestr(str2num(series_spec{2}(6:20))./10.^9)]),
        %
        %
        % plot complex
        position = [round(scrsz(3)/3.6*3),round(scrsz(4)*3/4)-18,...
            round(scrsz(3)/6),round(scrsz(4)/4-40)];
        h_complex = findobj('Position',position);
        if (isempty(h_complex)),
            h_complex = viewer_cmplx_fig;
        end,
        figure(h_complex), hold off,
        [H_atomic_real,H_atomic_imag]=pol2cart(angle(H_atomic(1:l))+pi/2,...
            20*log10(abs(H_atomic(1:l))/abs(H_atomic(F2(n)+1))*100000));
        compass(H_atomic_real(F2(n)+1),H_atomic_imag(F2(n)+1),'r'),
        hold on,zoom on,
        compass(H_atomic_real([F1(n)+1,F2(n)+1]),...
            H_atomic_imag([F1(n)+1,F2(n)+1]),'k'),
        if(L_TT(n)>10),
            compass(H_atomic_real(F_TT(n)+1),H_atomic_imag(F_TT(n)+1),'k'),
        end,
        color_order = [0 0 1;0 .75 0; 1 0 1; ...
            .75 .5 0; 1 .75 0;.75 .75 .75; 0 .75 .75; 0 0 0];
        for (o = 1:order),
            try
                h=compass(H_atomic_real(F2(n)+F_TT(n)*o+1),...
                    H_atomic_imag(F2(n)+F_TT(n)*o+1));
                set(h,'Color',color_order(o,:)),
                h=compass(H_atomic_real(F2(n)-F_TT(n)*o+1),...
                    H_atomic_imag(F2(n)-F_TT(n)*o+1));
                set(h,'Color',color_order(o,:)),
            catch
            end
        end,
        set(h_complex,'Position',position,...
            'Name',['DP_TT_cmplx: Nr.',num2str(n),', Ord. ',...
            num2str(order),', ',series_spec{1},' , ',...
            datestr(str2num(series_spec{2}(6:20))./10.^9)]),
        drawnow,


    case 'close_all'
        h = findobj('Type','figure','Tag','plot');
        if (isempty(h)),
            msgbox('No open plot !','Message: DP_TT Viewer','error'),
        else,
            close(h),
        end,

    case 'slider'
        val = get(findobj(gcf,'Tag','slider'),'Value');
        set(findobj(gcf,'Tag','order'),'String',num2str(round(val)));


    case 'close'
        h = findobj('Type','figure','Tag','plot');
        if (~isempty(h)),
            button = questdlg('Close all Plot figures?',...
                'DP_TT Viewer');
            if (strcmp(button,'Yes')),
                close(h),
            end,
            if (strcmp(button,'Cancel')),
                return,
            end,
        end,
        closereq,

    case 'comment'
        series_spec = get(gcbf,'User');
        comment = get(findobj(gcbf,'Tag','comment'),'String');
        save([OAE_PATH,'\Subjects\',series_spec{1},'\',series_spec{2}], ...
            'comment','-append'),
        dp_tt_browser('fill',series_spec{1}),


    case 'help'


end, % switch(action)


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
function DP_TT_viewer(action, patient, measurement)

global OAE_PATH POS_USED POS_TOGGLE
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
            'L1','l1','F2','L2','l2','F_TT','L_TT','l_TT','l_dp',...
            'comment','min_freq'),

        n = length(l1);
        if (~exist('l_dp')),
            l_dp=-99*ones(1,n);
        elseif size(l_dp,2) == 1
            l_dp = [l_dp' zeros(n,4)];
        end,
        set(findobj(h_fig,'Tag','comment'),'String',comment),
        if (length(l1)==0),
            return,
        end,
        if (~exist('L_TT_sl')),
            L_TT_sl = L_TT;
            L1_sl = L1;
            L2_sl = L2;
        end,
        string = ''; n_dp=0; n_dp2=0; n_sf=0;series_DP = [];series_DP2 = [];series_SF = [];
        for (n = 1:length(l1)),
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

            if L1(n)>= 0
                if F2(n)/F1(n) > 1.2
                    n_dp = n_dp+1;
                    series_DP(n_dp,:) = [n, floor(F_TT(n)*min_freq),l_TT(n),L_TT(n),floor(F1(n)*min_freq),...
                        l1(n),L1(n),floor(F2(n)*min_freq),l2(n),L2(n),l_dp(n,:)];
                else
                    n_dp2 = n_dp2+1;
                    series_DP2(n_dp2,:) = [n, floor(F_TT(n)*min_freq),l_TT(n),L_TT(n),floor(F1(n)*min_freq),...
                        l1(n),L1(n),floor(F2(n)*min_freq),l2(n),L2(n),l_dp(n,:)];
                end
            else
                n_sf = n_sf+1;
                series_SF(n_sf,:) = [n, floor(F_TT(n)*min_freq),l_TT(n),L_TT(n),floor(F1(n)*min_freq),...
                    l1(n),L1(n),floor(F2(n)*min_freq),l2(n),L2(n),l_dp(n,1),l_dp(n,2),...
                    l_dp(n,3),l_dp(n,4),l_dp(n,5),l_dp(n,6)];
            end

%             string{n}=sprintf(....
%                 '%2d%3d%1s%5d%1s%5d%1s%5.1f%5.1f%5.1f%5.1f%5.1f%5.1f',...
%                 n, floor(F_TT(n)*min_freq),str_TT,floor(F1(n)*min_freq),...
%                 str_L1,floor(F2(n)*min_freq),str_L2,l_dp(n,1),l_dp(n,2),...
%                 l_dp(n,3),l_dp(n,4),l_dp(n,5),l_dp(n,6));
            string{n}=sprintf(....
                '%2d%3d%1s%5d%1s%5d%1s%5.1f%5.1f%5.1f%5.1f%5.1f%5.1f',...
                n, floor(F_TT(n)*min_freq),str_TT,floor(F1(n)*min_freq),...
                str_L1,floor(F2(n)*min_freq),str_L2,l_dp(n,1),l_dp(n,2),...
                l_dp(n,3),l_dp(n,4),l_dp(n,5));
        end
        h = findobj(h_fig,'Tag','atomics');
        set(h,'String',string),
        set(h,'Value',n);
        format SHORTG

        if ~isempty(series_DP2)
            position = [680   260   560   210];
            h_series_DP = findobj('Position',position);
            if (isempty(h_series_DP))
                h_series_DP = figure;
            end
            figure(h_series_DP), clf,
            set(h_series_DP,'Position',position,...
                'Name',['DP_series: ',series_spec{1}(1:3),...
                ' , ',series_spec{2}(6:20)])
            hold off
            plot(l_dp(3:2:end,3),'bx-')
            hold on
            plot(l_dp(3:2:end,4),'rx-')
            plot(l_dp(3:2:end,5),'rx-')
            plot(l_dp(3:2:end,6),'bx-'),
            grid on, axis([0 75 -30 30]),title modSF
            %legend('-2','-1','+1','+2','Location','southeast')
        end

        if ~isempty(series_DP)
            position = [680   470   560   210];
            h_series_DP = findobj('Position',position);
            if (isempty(h_series_DP))
                h_series_DP = figure;
            end
            figure(h_series_DP), clf,
            set(h_series_DP,'Position',position,...
                'Name',['DP_series: ',series_spec{1}(1:3),...
                ' , ',series_spec{2}(6:20)])
            hold off
            plot(l_dp(2:2:end,1),'kx-', 'LineWidth', 2),
            hold on
            plot(l_dp(2:2:end,3),'bx-')
            plot(l_dp(2:2:end,4),'rx-')
            plot(l_dp(2:2:end,5),'rx-')
            plot(l_dp(2:2:end,6),'bx-')
            plot(l_dp(2:2:end,2),'ko:', 'LineWidth', 2),
            plot(l_dp(2:2:end,7),'bo:')
            plot(l_dp(2:2:end,8),'ro:')
            plot(l_dp(2:2:end,9),'ro:')
            plot(l_dp(2:2:end,10),'bo:')
            grid on, axis([0 75 -20 40]),title modDP, xlabel('no of 5s-presentions')
            %legend('2f1-f2','-2','-1','+1','+2','f2-f1','-2','-1','+1','+2','Location','southeast')
        end
        if ~isempty(series_SF)
            position = [680   678   560   210];
            h_series_SF = findobj('Position',position);
            if (isempty(h_series_SF))
                h_series_SF = figure;
            end
            figure(h_series_SF), clf,
            set(h_series_SF,'Position',position,...
                'Name',['SF_series: ',series_spec{1}(1:3),...
                ' , ',series_spec{2}(6:20)])
            plot(series_SF(:,13:16),'x-'), grid on
            %legend('-2','-1','+1','+2','Location','northeast')
        end



    case 'plot_replace',
        if(exist('patient')),
            h_fig = findobj('Name',['DP_TT_Viewer: ', patient,' , ',measurement(6:20)]);
        else,
            h_fig = gcbf;
        end,
        scrsz = get(0,'ScreenSize');
        n = get(findobj(h_fig,'Tag','atomics'),'Value');
        % get order of modulation
        order = str2num(get(findobj(h_fig,'Tag','order'),'String'));
        if (isempty(order)),
            msgbox('define the order of modulation',...
                'ERROR: DP_TT Viewer','error')
            return,
        end
        series_spec = get(h_fig,'User');
        load([OAE_PATH,'\Subjects\',series_spec{1},'\',series_spec{2}],...
            'min_freq', 'F1','F2','F_TT','L_TT','l_dp'),

        load([OAE_PATH,'\Subjects\',series_spec{1},'\', series_spec{2}], ...
            ['avg', num2str(n)], 'h_mic'),
        eval(['H_atomic = fft(avg',num2str(n), ')./h_mic;']);
        % Use for old measurments of eg. Manuela Berger und Jana Wolf
        % load([OAE_PATH,'\OAE\Subjects\',series_spec{1},'\', series_spec{2}],'gains'),
        % eval(['H_atomic = (fft(-avg',num2str(n),'./10^(gains(n)/20))./h_mic)'';']);
        % AND REPLACE line 392
        %    course = abs(ifft(H_atomic .* mod_mask));
        %    course = abs(ifft((H_atomic .* mod_mask)'));
        % i.e. IFFT gives reversed time course when applying to transposed matrix!
        % END: Use for old measurments of eg. Manuela Berger und Jana Wolf

        l = length(H_atomic);

        mod_mask(l,1) = 0;
        noise_mask1(l,1) = 0;
        noise_mask2(l,1) = 0;
        f(1) = 2*F1(n)-F2(n)+1;
        f(2) = F2(n)-2*F1(n)+1;
        f(3) = F2(n)-F1(n)+1;
        f(4) = F2(n)+F1(n)+1;
        f(5) = F2(n)+2*F1(n)+1;
        for q=1:5
            if f(q)>0 & f(q)< l
                mod_mask(f(q)) = 1;
            end
        end
        noise_mask1(l,1) = 0;
        noise_mask2(l,1) = 0;

        % plot spectrum
        position = [0,round(scrsz(4)*3/4)-18,round(scrsz(3)/3.6),...
            round(scrsz(4)/4-40)];
        h_spectrum = findobj('Position',position);
        if (isempty(h_spectrum)),
            h_spectrum = viewer_spctr_fig;
        end,
        figure(h_spectrum), clf,
        % plot spectrum
        x = [0:min_freq:l*min_freq];
        hold on,

        stem_60(x(1:end-1),20*log10(abs(H_atomic(1:length(x)-1))),'kx'),
        % stem_60(x(1:F_TT(n):length(x)),20*log10(abs(H_atomic(1:F_TT(n):length(x)))),'.'),

        stem_60(x([F1(n)+1,F2(n)+1]),20*log10(abs(H_atomic...
            ([F1(n)+1,F2(n)+1]))),'filled','bo'),

        stem_60(x(mod_mask>0),20*log10(abs(H_atomic(mod_mask>0))),'filled','ro'),

        f1_harmonics_mask = [];
        for q = 1: 1400/F1(n)
            f1_harmonics_mask = [f1_harmonics_mask q*F1(n)];
        end

        stem_60(x(f1_harmonics_mask+1),20*log10(abs(H_atomic...
            (f1_harmonics_mask+1))),'bx'),

        axis([0 24000 -40 115]),grid on, zoom on,
        set(h_spectrum,'Position',position,...
            'Name',['DP_TT_spctr: Nr.',num2str(n),', ',series_spec{1}(1:3),...
            ' , ',series_spec{2}(6:20)]),

        if ~isinf(F_TT(n)) & L_TT(n)> 70
            % plot magnitude
            position = [round(scrsz(3)/3.6),round(scrsz(4)*3/4)-18,...
                round(scrsz(3)/4),round(scrsz(4)/4-40)];
            h_magnitude = findobj('Position',position);
            if (isempty(h_magnitude)),
                h_magnitude = viewer_magn_fig;
            end,
            figure(h_magnitude), hold off,
            x = linspace(0,2/F_TT(n)/min_freq*1000,round(2*l/F_TT(n)));
            % plot time course of DP magintude
            hold off
            if F1(n)*min_freq > 500
                mod_mask = zeros(l,1);
                mod_mask([...
                    2*F1(n)-F2(n)+1, ...
                    2*F1(n)-F2(n)+1-F_TT(n), ...
                    2*F1(n)-F2(n)+1+F_TT(n), ...
                    %                 2*F1(n)-F2(n)+1+2*F_TT(n), ...
                    %                 2*F1(n)-F2(n)+1+F_TT(n), ...
                    ]) = 1;
                course = abs(ifft((H_atomic .* mod_mask)));
                plot(x,20*log10(l * course(1:round(2*l/F_TT(n))))),
                hold on;
                mod_mask = zeros(l,1);
                mod_mask([...
                    F2(n)-F1(n)+1, ...
                    F2(n)-F1(n)+1-F_TT(n), ...
                    F2(n)-F1(n)+1+F_TT(n), ...
                    %                 F2(n)-F1(n)+1+2*F_TT(n), ...
                    %                 F2(n)-F1(n)+1+F_TT(n), ...
                    ]) = 1;
                course = abs(ifft((H_atomic .* mod_mask)));
                plot(x,20*log10(l * course(1:round(2*l/F_TT(n)))),'b--'),
                legend({'2f2-f1' 'f2-f1'})

                if(L_TT(n)>70)
                    H_f_TT(l) = 0;
                    H_f_TT(F_TT(n)+1) = H_atomic(F_TT(n)+1);
                    course_TT = real(ifft(H_f_TT));
                    course_TT = course_TT/max(course_TT);
                    plot(x,20*course_TT(1:round(2*l/F_TT(n))),'k:'),
                end
                axis([0 2/F_TT(n)/min_freq*1000 -20 20]), grid on, zoom on,
            else
                f(1) = F2(n)-2*F1(n)+1;
                f(2) = F2(n)-F1(n)+1;
                f(3) = F2(n)+F1(n)+1;
                f(4) = F2(n)+2*F1(n)+1;
                for q=1:4
                    mod_mask = zeros(l,1);
                    mod_mask([...
                        f(q),...
                        f(q)-F_TT(n), ...
                        f(q)+F_TT(n), ...
                        ]) = 1;

                    course = abs(ifft((H_atomic .* mod_mask)));
                    plot(x,20*log10(l * course(1:round(2*l/F_TT(n))))),
                    hold on;
                end
                legend({'-2' '-1' '+1' '+2'})
                if(L_TT(n)>70)
                    H_f_TT(l) = 0;
                    H_f_TT(F_TT(n)+1) = H_atomic(F_TT(n)+1);
                    course_TT = real(ifft(H_f_TT));
                    course_TT = course_TT/max(course_TT);
                    plot(x,10*course_TT(1:round(2*l/F_TT(n)))+5,'k:'),
                end
                axis([0 2/F_TT(n)/min_freq*1000 -10 20]), grid on, zoom on,
            end

            %         MindB = min(20*log10(l * course(1:round(2*l/F_TT(n)))));
            %         plot(x, MindB*ones(1,length(((1:round(2*l/F_TT(n)))))),'Color','red','LineStyle','--'),
            %annotation('textbox', [0.2 0.1 0.3 0.3],'String',num2str(MindB),'FitBoxToText','on'),
            %legend(num2str(MindB),'Location','northwest')

            %     %plot noise  MODIFIED98-11-12H displ. of Noise in LINES 381-385
            %     if henselnoise, cour1=ifft(H_atomic.*noise_mask1);cour2=ifft(H_atomic.*noise_mask2);course=abs(cour1+cour2)/2;
            %     else course1 = abs(ifft(H_atomic .* noise_mask1));
            %         course2 = abs(ifft(H_atomic .* noise_mask2));
            %         course = (course1 + course2)/2;
            %     end,
            %     plot(x,20*log10(l * course(1:round(2*l/F_TT(n)))),':'),

            %     MindB = min(20*log10(l * course(1:round(2*l/F_TT(n)))));
            %     plot(x, MindB*ones(1,length(((1:round(2*l/F_TT(n)))))),'Color','red','LineStyle','--');
            %
            % plot time course of  the masker :387:
            %  time course of masker (scaled to 1)
            set(h_magnitude,'Position',position,...
                'Name',['DP_TT_mag: Nr.',num2str(n),', ',series_spec{1}(1:3),' , ',...
                series_spec{2}(6:20)]),
        end
        %
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
        % 	% plot complex
        %     position = [round(scrsz(3)/3.6*3),round(scrsz(4)*3/4)-18,...
        %             round(scrsz(3)/6),round(scrsz(4)/4-40)];
        %     h_complex = findobj('Position',position);
        %     if (isempty(h_complex)),
        %         h_complex = viewer_cmplx_fig;
        % 	end,
        %     figure(h_complex), hold off,
        %     [H_atomic_real,H_atomic_imag]=pol2cart(angle(H_atomic(1:l))+pi/2,...
        %         20*log10(abs(H_atomic(1:l))/abs(H_atomic(f_DP+1))*100000));
        %     compass(H_atomic_real(f_DP+1),H_atomic_imag(f_DP+1),'r'),
        %     hold on,zoom on,
        %     compass(H_atomic_real([F1(n)+1,F2(n)+1]),...
        %         H_atomic_imag([F1(n)+1,F2(n)+1]),'k'),
        %     if(L_TT(n)>50),
        %         compass(H_atomic_real(F_TT(n)+1),H_atomic_imag(F_TT(n)+1),'k'),
        %     end,
        %     color_order = [0 0 1;0 .75 0; 1 0 1; ...
        %            .75 .5 0; 1 .75 0;.75 .75 .75; 0 .75 .75; 0 0 0];
        %     for (o = 1:order),
        %         h=compass(H_atomic_real(f_DP+F_TT(n)*o+1),...
        %             H_atomic_imag(f_DP+F_TT(n)*o+1));
        %         set(h,'Color',color_order(o,:)),
        %         h=compass(H_atomic_real(f_DP-F_TT(n)*o+1),...
        %             H_atomic_imag(f_DP-F_TT(n)*o+1));
        %         set(h,'Color',color_order(o,:)),
        % 	end,
        % 	set(h_complex,'Position',position,...
        % 		'Name',['DP_TT_cmplx: Nr.',num2str(n),', Ord. ',...
        % 		num2str(order),', ',series_spec{1}(1:3),' , ',...
        % 		datestr(str2num(series_spec{2}(6:20))./10.^9)]),
        % 	drawnow,


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
        save([OAE_PATH,'\OAE\Subjects\',series_spec{1},'\',series_spec{2}], ...
            'comment','-append'),
        dp_tt_browser('fill',series_spec{1}),


    case 'help'


end, % switch(action)


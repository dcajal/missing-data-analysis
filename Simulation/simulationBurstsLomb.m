%% Lomb's method results with burst of missed beats

clear
addpath('lib','Simulation/database');
figurePresets
rng('default')


% Spectral analysis
windowSeconds = 60;
overlapSeconds = 30;
nfft = 2^10;

% subject = {'01M'};
subject = {'01M' '02V' '04V' '05V' '06M' '07M' '11M' '13V' '17V'}; % Respiración en banda de HF
% subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '10V' '11M'...
%     '12V' '13V' '15V' '16V' '17V'};
burstDuration = [0 5 10 15 20];

fillGaps = 'removeOutliers'; % 'removeOutliers' 'iterative' 'iterativeNonLinear'
nRealizations = 10;
displayCounter = 1;
plotflag = false;

resultsSupineLomb = cell(length(subject)*nRealizations,length(burstDuration));
resultsTiltLomb = cell(length(subject)*nRealizations,length(burstDuration));
for kk = 1:length(subject)
    for ll = 1:nRealizations
        for jj = 1:length(burstDuration)
            fprintf('Computing realization %i of Subject %s with %i burst duration (%i/%i)...',...
                ll,subject{kk},burstDuration(jj),displayCounter,...
                length(subject)*nRealizations*length(burstDuration));
            [SupineECG, TiltECG] = loadSimulationSubject(pwd,subject{kk});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Pulse deletion
            if burstDuration(jj) > 0
                SupineECG.tk(find(SupineECG.tk<=30+(60-burstDuration(jj))/(nRealizations-1)*(ll-1),1,'Last'):...
                    find(SupineECG.tk>=30+(60-burstDuration(jj))/(nRealizations-1)*(ll-1)+burstDuration(jj),1)) = [];
                TiltECG.tk(find(TiltECG.tk<=30+(60-burstDuration(jj))/(nRealizations-1)*(ll-1),1,'Last'):...
                    find(TiltECG.tk>=30+(60-burstDuration(jj))/(nRealizations-1)*(ll-1)+burstDuration(jj),1)) = [];

                switch fillGaps
                    case 'iterativeNonLinear'
                        SupineECG.tn = gapcorrectorNonLinear(SupineECG.tk); SupineECG.ids = [];
                        TiltECG.tn = gapcorrectorNonLinear(TiltECG.tk); TiltECG.ids = [];
                    case 'iterative'
                        SupineECG.tn = gapcorrector(SupineECG.tk);
                        TiltECG.tn = gapcorrector(TiltECG.tk);
                    case 'incidences'
                        [~,~,~,SupineECG.tn] = incidences(SupineECG.tk);
                        [~,~,~,TiltECG.tn] = incidences(TiltECG.tk);
                    case 'removeOutliers'
                        SupineECG.tn = SupineECG.tk;
                        TiltECG.tn = TiltECG.tk;
                end
            else
                SupineECG.tn = SupineECG.tk;
                TiltECG.tn = TiltECG.tk;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Lomb periodogram
            % Supine
            SupineLomb.f = linspace(0,0.4,nfft);
            SupineLomb.dtk = diff(SupineECG.tn);
            SupineLomb.tk = SupineECG.tn(2:end);

            % Delete peaks
            if strcmp(fillGaps,'removeOutliers')
                threshold = computeThreshold(SupineLomb.dtk);
                SupineLomb.tk(SupineLomb.dtk>threshold) = [];
                SupineLomb.dtk(SupineLomb.dtk>threshold) = [];
            end
            
            printAux = SupineLomb.dtk; 
            SupineLomb.hr = 1000./SupineLomb.dtk;
            SupineLomb.hr = detrend(SupineLomb.hr);
            clear threshold

            SupineLomb.pxx = plombav(SupineLomb.tk,SupineLomb.hr,SupineLomb.f,windowSeconds,overlapSeconds);

             [bb, aa] = butter(15, 0.04*2, 'high');
            h = freqz(bb,aa,nfft);
            SupineLomb.pxx = SupineLomb.pxx.*abs(h);
            SupineLomb.pxx = SupineLomb.pxx;
            clear h bb aa


            % Tilt
            TiltLomb.f = linspace(0,0.4,nfft);
            TiltLomb.dtk = diff(TiltECG.tn);
            TiltLomb.tk = TiltECG.tn(2:end);

            % Delete peaks
            if strcmp(fillGaps,'removeOutliers')
                threshold = computeThreshold(TiltLomb.dtk);
                TiltLomb.tk(TiltLomb.dtk>threshold) = [];
                TiltLomb.dtk(TiltLomb.dtk>threshold) = [];
            end
            TiltLomb.hr = 1000./TiltLomb.dtk;
            TiltLomb.hr = detrend(TiltLomb.hr);
            clear threshold

            TiltLomb.pxx = plombav(TiltLomb.tk,TiltLomb.hr,TiltLomb.f,windowSeconds,overlapSeconds);

            [bb, aa] = butter(15, 0.04*2, 'high');
            h = freqz(bb,aa,nfft);
            TiltLomb.pxx = TiltLomb.pxx.*abs(h);
            TiltLomb.pxx = TiltLomb.pxx;
            clear h bb aa           

    %       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plots

            if plotflag
                if burstDuration(jj)==0
                    figure('DefaultAxesFontSize',14,'units','normalized','outerposition',[0 0 1 1]); 
                    subplot(211)
                    stem(SupineLomb.tk,printAux,'LineWidth',2); hold on
                    ylabel('NN [s]');
                    xlabel('Time [s]');
                    title(subject{kk})
                    subplot(212)
                    plot(SupineLomb.f,SupineLomb.pxx,'LineWidth',2); hold on
                    xlabel('Frequency [Hz]')
                    ylabel('PSD_{Lomb} [ms^2]','interpreter','tex')
                    axis tight
                else
                    subplot(211)
                    stem(SupineLomb.tk,printAux,'LineWidth',1); hold on
                    subplot(212)
                    plot(SupineLomb.f,SupineLomb.pxx,'LineWidth',1);
                end
                if burstDuration(jj)==20
                    legend('0','5','10','15','20');
        %             set(gcf,'position',[0,0,1000,300])
        %             set(gcf, 'Position', get(0, 'Screensize'));
                    pause
                end
            end
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            resultsSupineLomb{(ll-1)*length(subject)+kk,jj} = freqind(SupineLomb.pxx,SupineLomb.f);
            resultsTiltLomb{(ll-1)*length(subject)+kk,jj} = freqind(TiltLomb.pxx,TiltLomb.f);
            
            fprintf('Done\n');  
            displayCounter = displayCounter+1;
        end
    end
end

clear SupineECG SupineLomb TiltECG TiltLomb jj kk fs nfft overlapSeconds windowSeconds

%% Degradation results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Degradation results (error): Frequency indexes (Lomb), error bursts, using %s method\n',fillGaps);
disp('Measure          Burst duration (seconds)');
fprintf('               '); fprintf('%i                    ',burstDuration(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

error = computeError(resultsSupineLomb, resultsTiltLomb, burstDuration, 'LF');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('P_{LF} [mHz^2]','interpreter','tex')
fprintf('P_LF           '); fprintf('%.2f (%.2f -- %.2f)   & ',error); fprintf('\n')

error = computeError(resultsSupineLomb, resultsTiltLomb, burstDuration, 'HF');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('P_{HF} [mHz^2]','interpreter','tex')
fprintf('P_HF            '); fprintf('%.2f (%.2f -- %.2f)   & ',error); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')

clear error

%% Sympathovagal balance results
 
fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('p-values (supine/tilt groups): Frequency indexes (Lomb), error bursts, using %s method\n',fillGaps);
disp('Measure        Burst duration (seconds)');
fprintf('               '); fprintf('%i          ',burstDuration(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

% figure(1)
% significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, burstDuration, 'LF');
% ylabel('P_{LF} [mHz^2]','interpreter','tex')
% fprintf('P_LF           '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

figure(2)
significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, burstDuration, 'HF');
ylabel('P_{HF} [mHz^2]','interpreter','tex')
fprintf('P_HF           '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

figure(3)
significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, burstDuration, 'LFn');
ylabel('P_{LFn}','interpreter','tex')
fprintf('P_LFn          '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

figure(4)
significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, burstDuration, 'LFHF');
ylabel('P_{LF}/P_{HF}','interpreter','tex')
fprintf('P_LF/P_HF      '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

clear significance

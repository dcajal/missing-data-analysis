%% Lomb's method results with burst of missed beats

clear
addpath('lib','Simulation/database');
figurePresets

% Spectral analysis
windowSeconds = 60;
overlapSeconds = 30;
nfft = 2^10;

% subject = {'07M'};
subject = {'01M' '02V' '04V' '05V' '06M' '07M' '11M' '13V' '17V'}; % Respiración en banda de HF
% subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '11M'...
%     '12V' '13V' '15V' '16V' '17V'};
burstDuration = [0 10 20 30 40 50 60];
fillGaps = 'iterativeNonLinear'; % 'removeOutliers' 'incidences' 'iterative' 'iterativeNonLinear'

resultsSupineLomb = cell(length(subject),length(burstDuration));
resultsTiltLomb = cell(length(subject),length(burstDuration));
for kk = 1:length(subject)
    for jj = 1:length(burstDuration)
        [SupineECG, TiltECG] = loadSimulationSubject(pwd,subject{kk});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Pulse deletion
        if burstDuration(jj) > 0
            rng('default')
    
            SupineECG.tk(find(SupineECG.tk<=60-burstDuration(jj)/2,1,'Last'):...
                find(SupineECG.tk>=60+burstDuration(jj)/2,1)) = [];
            TiltECG.tk(find(TiltECG.tk<=60-burstDuration(jj)/2,1,'Last'):...
                find(TiltECG.tk>=60+burstDuration(jj)/2,1)) = [];

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
        SupineLomb.dtk = 1000.*SupineLomb.dtk;
        SupineLomb.dtk = detrend(SupineLomb.dtk);
        clear threshold

        SupineLomb.pxx = plombav(SupineLomb.tk,SupineLomb.dtk,SupineLomb.f,windowSeconds,overlapSeconds);

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
        TiltLomb.dtk = 1000.*TiltLomb.dtk;
        TiltLomb.dtk = detrend(TiltLomb.dtk);
        clear threshold

        TiltLomb.pxx = plombav(TiltLomb.tk,TiltLomb.dtk,TiltLomb.f,windowSeconds,overlapSeconds);
        
        [bb, aa] = butter(15, 0.04*2, 'high');
        h = freqz(bb,aa,nfft);
        TiltLomb.pxx = TiltLomb.pxx.*abs(h);
        TiltLomb.pxx = TiltLomb.pxx;
        clear h bb aa           
        
%         if jj==1
%             figure('DefaultAxesFontSize',14); hold on;
%             plot(TiltLomb.f,TiltLomb.pxx,'LineWidth',2); hold on
%             title(subject{kk})
%         else
%             plot(TiltLomb.f,TiltLomb.pxx,'LineWidth',1); hold on
%         end
%         xlabel('Frequency [Hz]')
%         ylabel('PSD_{Lomb} [ms^2/Hz]','interpreter','tex')
%         axis tight
%         set(gcf,'position',[0,0,1000,300])
%         if jj==length(burstDuration)
%             legend('0s','10s','20s','30s','40s','50s','60s');
%         end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        resultsSupineLomb{kk,jj} = freqind(SupineLomb.pxx,SupineLomb.f);
        resultsTiltLomb{kk,jj} = freqind(TiltLomb.pxx,TiltLomb.f);
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
% ylabel('P_{LF} [ms^2]','interpreter','tex')
fprintf('P_LF           '); fprintf('%.2f (%.2f-%.2f)   ',error); fprintf('\n')


error = computeError(resultsSupineLomb, resultsTiltLomb, burstDuration, 'HF');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('P_{HF} [ms^2]','interpreter','tex')
fprintf('P_HF            '); fprintf('%.2f (%.2f-%.2f)   ',error); fprintf('\n')


error = computeError(resultsSupineLomb, resultsTiltLomb, burstDuration, 'LFn');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('P_{LFn}','interpreter','tex')
fprintf('P_LFn          '); fprintf('%.2f (%.2f-%.2f)       ',error); fprintf('\n')


error = computeError(resultsSupineLomb, resultsTiltLomb, burstDuration, 'LFHF');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('P_{LF}/P_{HF}','interpreter','tex')
fprintf('P_LF/P_HF      '); fprintf('%.2f (%.2f-%.2f)       ',error); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')

clear error

%% Sympathovagal balance results
 
fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('p-values (supine/tilt groups): Frequency indexes (Lomb), error bursts, using %s method\n',fillGaps);
disp('Measure        Burst duration (seconds)');
fprintf('               '); fprintf('%i          ',burstDuration(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, burstDuration, 'LF');
% ylabel('P_{LF} [ms^2]','interpreter','tex')
fprintf('P_LF           '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, burstDuration, 'HF');
% ylabel('P_{HF} [ms^2]','interpreter','tex')
fprintf('P_HF           '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, burstDuration, 'LFn');
% ylabel('P_{LFn}','interpreter','tex')
fprintf('P_LFn          '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, burstDuration, 'LFHF');
% ylabel('P_{LF}/P_{HF}','interpreter','tex')
fprintf('P_LF/P_HF      '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

clear significance

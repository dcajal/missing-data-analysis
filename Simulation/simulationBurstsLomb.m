%% Lomb's method results with burst of missed beats

clear
addpath('lib','Simulation/database');
figurePresets

% Spectral analysis
windowSeconds = 60;
overlapSeconds = 30;
nfft = 2^10;

% subject = {'07M'};
% subject = {'01M' '02V' '04V' '05V' '06M' '07M' '11M' '13V' '17V'}; % Respiración en banda de HF
subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '11M'...
    '12V' '13V' '15V' '16V' '17V'};
burstDuration = [0 10 20 30 40 50 60];
fillGaps = 'iterativeNonLinear'; % 'none' 'incidences' 'iterative' 'iterativeNonLinear'

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
                case 'none'
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
        if strcmp(fillGaps,'none')
            threshold = computeThreshold(SupineLomb.dtk);
            SupineLomb.tk(SupineLomb.dtk>threshold) = [];
            SupineLomb.dtk(SupineLomb.dtk>threshold) = [];
        end
        SupineLomb.dtk = 1000.*SupineLomb.dtk;
        SupineLomb.dtk = detrend(SupineLomb.dtk);
        clear threshold

        SupineLomb.pxx = plombav(SupineLomb.tk,SupineLomb.dtk,SupineLomb.f,windowSeconds,overlapSeconds);
%         SupineLomb.pxx = plomb(SupineLomb.dtk,SupineLomb.tk,SupineLomb.f);

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
        if strcmp(fillGaps,'none')
            threshold = computeThreshold(TiltLomb.dtk);
            TiltLomb.tk(TiltLomb.dtk>threshold) = [];
            TiltLomb.dtk(TiltLomb.dtk>threshold) = [];
        end
        TiltLomb.dtk = 1000.*TiltLomb.dtk;
        TiltLomb.dtk = detrend(TiltLomb.dtk);
        clear threshold

        TiltLomb.pxx = plombav(TiltLomb.tk,TiltLomb.dtk,TiltLomb.f,windowSeconds,overlapSeconds);
%         TiltLomb.pxx = plomb(TiltLomb.dtk,TiltLomb.tk,TiltLomb.f);

        
        [bb, aa] = butter(15, 0.04*2, 'high');
        h = freqz(bb,aa,nfft);
        TiltLomb.pxx = TiltLomb.pxx.*abs(h);
        TiltLomb.pxx = TiltLomb.pxx;
        clear h bb aa           
        
        if jj==1
            figure('DefaultAxesFontSize',14); hold on;
            plot(TiltLomb.f,TiltLomb.pxx,'LineWidth',2); hold on
            title(subject{kk})
        else
            plot(TiltLomb.f,TiltLomb.pxx,'LineWidth',1); hold on
        end
        xlabel('Frequency [Hz]')
        ylabel('PSD_{Lomb} [ms^2/Hz]','interpreter','tex')
        axis tight
        set(gcf,'position',[0,0,1000,300])
        if jj==length(burstDuration)
            legend('0s','10s','20s','30s','40s','50s','60s');
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % HF
        HFBand = [find(SupineLomb.f>=0.15,1) find(SupineLomb.f>=0.4,1)];
        SupineLomb.hf = trapz(SupineLomb.f(HFBand(1):HFBand(2)),SupineLomb.pxx(HFBand(1):HFBand(2)));
        TiltLomb.hf = trapz(TiltLomb.f(HFBand(1):HFBand(2)),TiltLomb.pxx(HFBand(1):HFBand(2)));
        
        % LF
        LFBand = [find(SupineLomb.f>=0.04,1) find(SupineLomb.f>=0.15,1)];
        SupineLomb.lf = trapz(SupineLomb.f(LFBand(1):LFBand(2)),SupineLomb.pxx(LFBand(1):LFBand(2)));
        TiltLomb.lf = trapz(TiltLomb.f(LFBand(1):LFBand(2)),TiltLomb.pxx(LFBand(1):LFBand(2)));
        
        % LFn
        SupineLomb.lfn = SupineLomb.lf/(SupineLomb.lf+SupineLomb.hf);
        TiltLomb.lfn = TiltLomb.lf/(TiltLomb.lf+TiltLomb.hf);
        
        % LF/HF
        SupineLomb.lfhf = SupineLomb.lf/SupineLomb.hf;
        TiltLomb.lfhf = TiltLomb.lf/TiltLomb.hf;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        % Save results
        resultsSupineLomb{kk,jj} = SupineLomb;
        resultsTiltLomb{kk,jj} = TiltLomb;
    end
end

clear SupineECG SupineLomb TiltECG TiltLomb jj kk fs nfft overlapSeconds windowSeconds HFBand LFBand

%% Degradation results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Degradation results (RMSE): Frequency indexes (Lomb), error bursts, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i          ',burstDuration(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')


RMSE = computeError(resultsSupineLomb, resultsTiltLomb, burstDuration, 'lf', true);
xlabel('Burst duration (s)','interpreter','tex')
ylabel('P_{LF} [ms^2]','interpreter','tex')
fprintf('P_LF (norm)    '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')


RMSE = computeError(resultsSupineLomb, resultsTiltLomb, burstDuration, 'hf', true);
xlabel('Burst duration (s)','interpreter','tex')
ylabel('P_{HF} [ms^2]','interpreter','tex')
fprintf('P_HF (norm)    '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')


RMSE = computeError(resultsSupineLomb, resultsTiltLomb, burstDuration, 'lfn');
xlabel('Burst duration (s)','interpreter','tex')
ylabel('P_{LFn}','interpreter','tex')
fprintf('P_LFn          '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')


RMSE = computeError(resultsSupineLomb, resultsTiltLomb, burstDuration, 'lfhf');
xlabel('Burst duration (s)','interpreter','tex')
ylabel('P_{LF}/P_{HF}','interpreter','tex')
fprintf('P_LF/P_HF      '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')

clear RMSE

%% Sympathovagal balance results
 
fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('p-values (supine/tilt groups): Frequency indexes (Welch), random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i          ',burstDuration); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, burstDuration, 'lf');
ylabel('P_{LF} [ms^2]','interpreter','tex')
fprintf('P_LF (norm)    '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, burstDuration, 'hf');
ylabel('P_{HF} [ms^2]','interpreter','tex')
fprintf('P_HF (norm)    '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, burstDuration, 'lfn');
ylabel('P_{LFn}','interpreter','tex')
fprintf('P_LFn          '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, burstDuration, 'lfhf');
ylabel('P_{LF}/P_{HF}','interpreter','tex')
fprintf('P_LF/P_HF      '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

clear significance

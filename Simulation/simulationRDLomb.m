%% Lomb's method results with random distributed errors

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
deletionProbability = 0:0.05:0.25;
fillGaps = 'iterativeNonLinear'; % 'none' 'incidences' 'iterative' 'iterativeNonLinear'

resultsSupineLomb = cell(length(subject),length(deletionProbability));
resultsTiltLomb = cell(length(subject),length(deletionProbability));
for kk = 1:length(subject)
    for jj = 1:length(deletionProbability)
        [SupineECG, TiltECG] = loadSimulationSubject(pwd,subject{kk});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Pulse deletion
        if deletionProbability(jj) > 0
            rng('default')
    
            SupineECG.tk(binornd(1,deletionProbability(jj),length(SupineECG.tk),1)>0.5) = [];
            TiltECG.tk(binornd(1,deletionProbability(jj),length(TiltECG.tk),1)>0.5) = [];

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

        [bb, aa] = butter(15, 0.04*2, 'high');
        h = freqz(bb,aa,nfft);
        TiltLomb.pxx = TiltLomb.pxx.*abs(h);
        TiltLomb.pxx = TiltLomb.pxx;
        clear h bb aa        
        
%         if deletionProbability(jj)==0
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
%         if deletionProbability(jj)==0.25
%             legend('0%','5%','10%','15%','20%','25%');
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
fprintf('Degradation results (RMSE): Frequency indexes (Lomb), random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i          ',100*deletionProbability(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')


RMSE = computeError(resultsSupineLomb, resultsTiltLomb, 100*deletionProbability, 'LF', true);
ylabel('P_{LF} [ms^2]','interpreter','tex')
fprintf('P_LF (norm)    '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')


RMSE = computeError(resultsSupineLomb, resultsTiltLomb, 100*deletionProbability, 'HF', true);
ylabel('P_{HF} [ms^2]','interpreter','tex')
fprintf('P_HF (norm)    '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')


RMSE = computeError(resultsSupineLomb, resultsTiltLomb, 100*deletionProbability, 'LFn');
ylabel('P_{LFn}','interpreter','tex')
fprintf('P_LFn          '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')


RMSE = computeError(resultsSupineLomb, resultsTiltLomb, 100*deletionProbability, 'LFHF');
ylabel('P_{LF}/P_{HF}','interpreter','tex')
fprintf('P_LF/P_HF      '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')

clear RMSE

%% Sympathovagal balance results
 
fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('p-values (supine/tilt groups): Frequency indexes (Lomb), random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i          ',deletionProbability(2:end)*100); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, 100*deletionProbability, 'LF');
ylabel('P_{LF} [ms^2]','interpreter','tex')
fprintf('P_LF (norm)    '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, 100*deletionProbability, 'HF');
ylabel('P_{HF} [ms^2]','interpreter','tex')
fprintf('P_HF (norm)    '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, 100*deletionProbability, 'LFn');
ylabel('P_{LFn}','interpreter','tex')
fprintf('P_LFn          '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, 100*deletionProbability, 'LFHF');
ylabel('P_{LF}/P_{HF}','interpreter','tex')
fprintf('P_LF/P_HF      '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

clear significance
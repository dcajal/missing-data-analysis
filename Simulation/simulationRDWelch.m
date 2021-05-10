%% Welch's method results with random distributed errors

clear
addpath('lib','Simulation/database');
figurePresets

% Spectral analysis
windowSeconds = 60;
overlapSeconds = 30;
nfft = 2^10;
fs = 4;

% subject = {'10V'};
subject = {'01M' '02V' '04V' '05V' '06M' '07M' '11M' '13V' '17V'}; % Respiración en banda de HF
% subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '10V' '11M'...
%     '12V' '13V' '15V' '16V' '17V'};
deletionProbability = 0:0.05:0.25;
fillGaps = 'iterativeNonLinear'; % 'ipfm' 'incidences' 'iterative' 'iterativeNonLinear'

resultsSupineIPFM = cell(length(subject),length(deletionProbability));
resultsTiltIPFM = cell(length(subject),length(deletionProbability));
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
                    SupineECG.tn = gapcorrector(SupineECG.tk); SupineECG.ids = [];
                    TiltECG.tn = gapcorrector(TiltECG.tk); TiltECG.ids = [];
                case 'incidences'
                    [~,~,~,SupineECG.tn] = incidences(SupineECG.tk); SupineECG.ids = [];
                    [~,~,~,TiltECG.tn] = incidences(TiltECG.tk); TiltECG.ids = [];
                case 'ipfm'
                    % tn is actually tm
                    [SupineECG.tn,SupineECG.ids] = incidences(SupineECG.tk);
                    [TiltECG.tn,TiltECG.ids] = incidences(TiltECG.tk);
            end
        else
            SupineECG.tn = SupineECG.tk; SupineECG.ids = [];
            TiltECG.tn = TiltECG.tk; TiltECG.ids = [];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % IPFM
        SupineIPFM.t = 0:1/fs:120;
        SupineIPFM.mt = ipfm(SupineECG.tn,SupineECG.ids,SupineIPFM.t);

        TiltIPFM.t = SupineIPFM.t;
        TiltIPFM.mt = ipfm(TiltECG.tn,TiltECG.ids,TiltIPFM.t);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Welch
        wdw = hamming(windowSeconds*fs);
        noverlap = overlapSeconds*fs;

        SupineIPFM.f = linspace(0,0.4,nfft);
        TiltIPFM.f = linspace(0,0.4,nfft);

        [bb, aa] = butter(4, 0.04*2/fs, 'high');
        auxSupine = filtfilt(bb, aa, 1000./SupineIPFM.mt(SupineIPFM.mt>0));
        auxTilt = filtfilt(bb, aa, 1000./TiltIPFM.mt(TiltIPFM.mt>0));

        SupineIPFM.pxx = pwelch(auxSupine,wdw,noverlap,SupineIPFM.f,fs);
        TiltIPFM.pxx = pwelch(auxTilt,wdw,noverlap,TiltIPFM.f,fs);
        clear wdw noverlap bb aa auxTilt auxSupine
        
%         if deletionProbability(jj)==0
%             figure('DefaultAxesFontSize',14); hold on;
%             plot(SupineIPFM.f,SupineIPFM.pxx,'LineWidth',2);
%             title(subject{kk})
%         else
%             plot(SupineIPFM.f,SupineIPFM.pxx,'LineWidth',1); hold on
%         end
%         xlabel('Frequency [Hz]')
%         ylabel('PSD_{Welch} [ms^2/Hz]','interpreter','tex')
%         axis tight
%         set(gcf,'position',[0,0,1000,300])
%         if deletionProbability(jj)==0.25
%             legend('0%','5%','10%','15%','20%','25%');
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        resultsSupineIPFM{kk,jj} = freqind(SupineIPFM.pxx,SupineIPFM.f);
        resultsTiltIPFM{kk,jj} = freqind(TiltIPFM.pxx,TiltIPFM.f);
             
    end
end

clear SupineECG SupineIPFM TiltECG TiltIPFM jj kk fs nfft overlapSeconds windowSeconds

%% Degradation results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Degradation results: Frequency indexes (Welch), random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i          ',100*deletionProbability(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

error = computeError(resultsSupineIPFM, resultsTiltIPFM, 100*deletionProbability, 'LF');
% ylabel('P_{LF} [ms^2]','interpreter','tex')
fprintf('P_LF           '); fprintf('%.2f (%.2f-%.2f)       ',error); fprintf('\n')


error = computeError(resultsSupineIPFM, resultsTiltIPFM, 100*deletionProbability, 'HF');
% ylabel('P_{HF} [ms^2]','interpreter','tex')
fprintf('P_HF           '); fprintf('%.2f (%.2f-%.2f)       ',error); fprintf('\n')


error = computeError(resultsSupineIPFM, resultsTiltIPFM, 100*deletionProbability, 'LFn');
% ylabel('P_{LFn}','interpreter','tex')
fprintf('P_LFn          '); fprintf('%.2f (%.2f-%.2f)       ',error); fprintf('\n')


error = computeError(resultsSupineIPFM, resultsTiltIPFM, 100*deletionProbability, 'LFHF');
% ylabel('P_{LF}/P_{HF}','interpreter','tex')
fprintf('P_LF/P_HF      '); fprintf('%.2f (%.2f-%.2f)       ',error); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')

clear error

%% Sympathovagal balance results
 
fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('p-values (supine/tilt groups): Frequency indexes (Welch), random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i          ',deletionProbability(2:end)*100); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

significance = twoGroupsDegradation(resultsSupineIPFM, resultsTiltIPFM, 100*deletionProbability, 'LF');
% ylabel('P_{LF} [ms^2]','interpreter','tex')
fprintf('P_LF           '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineIPFM, resultsTiltIPFM, 100*deletionProbability, 'HF');
% ylabel('P_{HF} [ms^2]','interpreter','tex')
fprintf('P_HF           '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineIPFM, resultsTiltIPFM, 100*deletionProbability, 'LFn');
% ylabel('P_{LFn}','interpreter','tex')
fprintf('P_LFn          '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineIPFM, resultsTiltIPFM, 100*deletionProbability, 'LFHF');
% ylabel('P_{LF}/P_{HF}','interpreter','tex')
fprintf('P_LF/P_HF      '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

clear significance

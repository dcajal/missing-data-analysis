%% Welch's method results with random distributed errors

clear
addpath('lib','database');
set(0,'defaultAxesFontName', 'Helvetica');
set(0,'defaultTextFontName', 'Helvetica');
set(0,'defaulttextInterpreter','tex')
set(0,'defaultLegendInterpreter','tex');
set(0,'defaultAxesFontSize', 20);
set(0,'defaultTextFontSize', 20);
set(0,'defaultFigureRenderer', 'painters')

% Signal duration
SIGNAL_DURATION = 120; % [s]

% Spectral analysis
WINDOW_SECONDS = 60; % Window size [s]
OVERLAP_SECONDS = 30; % Overlapping [s]
NFFT = 2^10;

subject = {'10V'};
% subject = {'01M' '02V' '04V' '05V' '06M' '07M' '11M' '13V' '17V'}; % Respiración en banda de HF
% subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '10V' '11M'...
%     '12V' '13V' '15V' '16V' '17V'};
deletionProbability = 0:0.05:0.25;
fillGaps = 'iterative'; % 'ipfm' 'incidences' 'iterative'
saveFigs = false;

resultsSupineIPFM = cell(length(subject),length(deletionProbability));
resultsTiltIPFM = cell(length(subject),length(deletionProbability));
for kk = 1:length(subject)
    for jj = 1:length(deletionProbability)
        load(strcat(pwd,'\database\',subject{kk},'.mat'));
        if SIGNAL_DURATION < 240
            tiltMarks(2) = tiltMarks(1) + SIGNAL_DURATION;
            tiltMarks(4) = tiltMarks(3) + SIGNAL_DURATION;
        end
        
        % Split Supine-Tilt
        SupineECG.signal = Ecg.signal(tiltMarks(1)*Ecg.samplerate:tiltMarks(2)*Ecg.samplerate);
        SupineECG.t = 0:1/Ecg.samplerate:(length(SupineECG.signal)-1)/Ecg.samplerate;
        SupineECG.tk = Ecg.qrs(find(Ecg.qrs>=tiltMarks(1),1):find(Ecg.qrs<=tiltMarks(2),1,'last'));
        SupineECG.tk = SupineECG.tk-tiltMarks(1);
        SupineECG.samplerate = Ecg.samplerate;

        TiltECG.signal = Ecg.signal(tiltMarks(3)*Ecg.samplerate:tiltMarks(4)*Ecg.samplerate);
        TiltECG.t = 0:1/Ecg.samplerate:(length(TiltECG.signal)-1)/Ecg.samplerate;
        TiltECG.tk = Ecg.qrs(find(Ecg.qrs>=tiltMarks(3),1):find(Ecg.qrs<=tiltMarks(4),1,'last'));
        TiltECG.tk = TiltECG.tk-tiltMarks(3);
        TiltECG.samplerate = Ecg.samplerate;

        clear tiltMarks Ecg Resp

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Pulse deletion
        if deletionProbability(jj) > 0
            rng('default')
    
            SupineECG.tk(binornd(1,deletionProbability(jj),length(SupineECG.tk),1)>0.5) = [];
            TiltECG.tk(binornd(1,deletionProbability(jj),length(TiltECG.tk),1)>0.5) = [];

            switch fillGaps
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
        SupineIPFM.t = SupineECG.t;
        SupineIPFM.mt = ipfm(SupineECG.tn,SupineECG.ids,SupineIPFM.t);

        TiltIPFM.t = SupineECG.t;
        TiltIPFM.mt = ipfm(TiltECG.tn,TiltECG.ids,TiltIPFM.t);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Welch
        wdw = hamming(WINDOW_SECONDS*SupineECG.samplerate);
        noverlap = OVERLAP_SECONDS*SupineECG.samplerate;

        SupineIPFM.f = linspace(0,0.4,NFFT);
        TiltIPFM.f = linspace(0,0.4,NFFT);

        [bb, aa] = butter(4, 0.04*2/SupineECG.samplerate, 'high');
        auxSupine = filtfilt(bb, aa, 1000./SupineIPFM.mt(SupineIPFM.mt>0));
        auxTilt = filtfilt(bb, aa, 1000./TiltIPFM.mt(TiltIPFM.mt>0));

        SupineIPFM.pxx = pwelch(auxSupine,wdw,noverlap,SupineIPFM.f,SupineECG.samplerate);
        TiltIPFM.pxx = pwelch(auxTilt,wdw,noverlap,TiltIPFM.f,TiltECG.samplerate);
        clear wdw noverlap bb aa auxTilt auxSupine
        
        if deletionProbability(jj)==0
            figure('DefaultAxesFontSize',14); hold on;
            plot(SupineIPFM.f,SupineIPFM.pxx,'LineWidth',2);
            title(subject{kk})
        else
            plot(SupineIPFM.f,SupineIPFM.pxx,'LineWidth',1); hold on
        end
        xlabel('Frequency [Hz]')
        ylabel('PSD_{Welch} [ms^2/Hz]')
        axis tight
%         ylim([0 11000])
        set(gcf,'position',[0,0,1000,300])
        if deletionProbability(jj)==0.25
            legend('0%','5%','10%','15%','20%','25%');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        % HF
        HFBand = [find(SupineIPFM.f>=0.15,1) find(SupineIPFM.f>=0.4,1)];
        SupineIPFM.hf = trapz(SupineIPFM.f(HFBand(1):HFBand(2)),SupineIPFM.pxx(HFBand(1):HFBand(2)));
        TiltIPFM.hf = trapz(TiltIPFM.f(HFBand(1):HFBand(2)),TiltIPFM.pxx(HFBand(1):HFBand(2)));
        
        % LF
        LFBand = [find(SupineIPFM.f>=0.04,1) find(SupineIPFM.f>=0.15,1)];
        SupineIPFM.lf = trapz(SupineIPFM.f(LFBand(1):LFBand(2)),SupineIPFM.pxx(LFBand(1):LFBand(2)));
        TiltIPFM.lf = trapz(TiltIPFM.f(LFBand(1):LFBand(2)),TiltIPFM.pxx(LFBand(1):LFBand(2)));
        
        % LFn
        SupineIPFM.lfn = SupineIPFM.lf/(SupineIPFM.lf+SupineIPFM.hf);
        TiltIPFM.lfn = TiltIPFM.lf/(TiltIPFM.lf+TiltIPFM.hf);
        
        % LF/HF
        SupineIPFM.lfhf = SupineIPFM.lf/SupineIPFM.hf;
        TiltIPFM.lfhf = TiltIPFM.lf/TiltIPFM.hf;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Save results
        resultsSupineIPFM{kk,jj} = SupineIPFM;
        resultsTiltIPFM{kk,jj} = TiltIPFM;
             
    end
end

%% Degradation results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Degradation results (RMSE): Frequency indexes (Welch), random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i          ',deletionProbability(2:end)*100); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineIPFM{kk,jj}.lf; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltIPFM{kk-length(subject),jj}.lf; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('P_{LF} [ms^2]','interpreter','tex')
% RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)); %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_welch_' fillGaps '_lf']),'epsc'); end %#ok<*UNRCH>
fprintf('P_LF (norm)    '); fprintf('%.3f       ',nRMSE(2:end)); fprintf('\n')


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineIPFM{kk,jj}.hf; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltIPFM{kk-length(subject),jj}.hf; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('P_{HF} [ms^2]','interpreter','tex')
% RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)); %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_welch_' fillGaps '_hf']),'epsc'); end
fprintf('P_HF (norm)    '); fprintf('%.3f       ',nRMSE(2:end)); fprintf('\n')


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineIPFM{kk,jj}.lfn; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltIPFM{kk-length(subject),jj}.lfn; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('P_{LFn}','interpreter','tex')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))); % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_welch_' fillGaps '_lfn']),'epsc'); end
fprintf('P_LFn          '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineIPFM{kk,jj}.lfhf; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltIPFM{kk-length(subject),jj}.lfhf; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('P_{LF}/P_{HF}','interpreter','tex')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))); % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_welch_' fillGaps '_lfhf']),'epsc'); end
fprintf('P_LF/P_HF      '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')


%% Sympathovagal balance results

% aux1 = [];
% aux2 = [];
% significance = [];
% for jj=1:length(deletionProbability)
%     for kk=1:length(subject)
%         aux1(kk,jj) = resultsSupineIPFM{kk,jj}.lf;
%         aux2(kk,jj) = resultsTiltIPFM{kk,jj}.lf;
%     end
%     significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
% end
% fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
% ylabel('P_{LF} [ms^2]')
% if saveFigs, saveas(gcf,strcat(['images\ans_rd_welch_' fillGaps '_lf']),'epsc'); end
% 
% 
% aux1 = [];
% aux2 = [];
% significance = [];
% for jj=1:length(deletionProbability)
%     for kk=1:length(subject)
%         aux1(kk,jj) = resultsSupineIPFM{kk,jj}.hf;
%         aux2(kk,jj) = resultsTiltIPFM{kk,jj}.hf;
%     end
%     significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
% end
% fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
% ylabel('P_{HF} [ms^2]')
% if saveFigs, saveas(gcf,strcat(['images\ans_rd_welch_' fillGaps '_hf']),'epsc'); end
% 
% 
% aux1 = [];
% aux2 = [];
% significance = [];
% for jj=1:length(deletionProbability)
%     for kk=1:length(subject)
%         aux1(kk,jj) = resultsSupineIPFM{kk,jj}.lfn;
%         aux2(kk,jj) = resultsTiltIPFM{kk,jj}.lfn;
%     end
%     significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
% end
% fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
% ylabel('P_{LFn}')
% if saveFigs, saveas(gcf,strcat(['images\ans_rd_welch_' fillGaps '_lfn']),'epsc'); end
% 
% 
% aux1 = [];
% aux2 = [];
% significance = [];
% for jj=1:length(deletionProbability)
%     for kk=1:length(subject)
%         aux1(kk,jj) = resultsSupineIPFM{kk,jj}.lfhf; %#ok<*SAGROW>
%         aux2(kk,jj) = resultsTiltIPFM{kk,jj}.lfhf;
%     end
%     significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
% end
% fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
% ylabel('P_{LF}/P_{HF}')
% if saveFigs, saveas(gcf,strcat(['images\ans_rd_welch_' fillGaps '_lfhf']),'epsc'); end


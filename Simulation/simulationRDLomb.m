%% Lomb's method results with random distributed errors

clear

% Signal duration
SIGNAL_DURATION = 120; % [s]

% Spectral analysis
WINDOW_SECONDS = 60; % Window size [s]
OVERLAP_SECONDS = 30; % Overlapping [s]
NFFT = 2^10;

% subject = {'07M'};
% subject = {'01M' '02V' '04V' '05V' '06M' '07M' '11M' '13V' '17V'}; % Respiración en banda de HF
subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '11M'...
    '12V' '13V' '15V' '16V' '17V'};
deletionProbability = [0 0.05 0.1 0.15 0.2 0.25];
fillGaps = 'incidences'; % 'none' 'incidences' 'iterative'
saveFigs = false;

resultsSupineLomb = cell(length(subject),length(deletionProbability));
resultsTiltLomb = cell(length(subject),length(deletionProbability));
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
                    SupineECG.tn = gapcorrectorNonLinear(SupineECG.tk);
                    TiltECG.tn = gapcorrectorNonLinear(TiltECG.tk);
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
        SupineLomb.f = linspace(0,0.4,NFFT);
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

        SupineLomb.pxx = plombav(SupineLomb.tk,SupineLomb.dtk,SupineLomb.f,WINDOW_SECONDS,OVERLAP_SECONDS);
%         SupineLomb.pxx = plomb(SupineLomb.dtk,SupineLomb.tk,SupineLomb.f);

        [bb, aa] = butter(15, 0.04*2, 'high');
        h = freqz(bb,aa,NFFT);
        SupineLomb.pxx = SupineLomb.pxx.*abs(h);
        SupineLomb.pxx = SupineLomb.pxx;
        clear h bb aa

        % Tilt
        TiltLomb.f = linspace(0,0.4,NFFT);
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

        TiltLomb.pxx = plombav(TiltLomb.tk,TiltLomb.dtk,TiltLomb.f,WINDOW_SECONDS,OVERLAP_SECONDS);
%         TiltLomb.pxx = plomb(TiltLomb.dtk,TiltLomb.tk,TiltLomb.f);

        
        [bb, aa] = butter(15, 0.04*2, 'high');
        h = freqz(bb,aa,NFFT);
        TiltLomb.pxx = TiltLomb.pxx.*abs(h);
        TiltLomb.pxx = TiltLomb.pxx;
        clear h bb aa        
        
%         if deletionProbability(jj)==0
%             figure('DefaultAxesFontSize',14); hold on;
%             plot(TiltLomb.f,TiltLomb.pxx,'LineWidth',2); hold on
% %             title(subject{kk})
%         else
%             plot(TiltLomb.f,TiltLomb.pxx,'LineWidth',1); hold on
%         end
%         xlabel('Frequency [Hz]')
%         ylabel('PSD_{Lomb} [ms^2/Hz]')
%         axis tight
%         ylim([0 12000])
%         set(gcf,'position',[0,0,1000,300])
%         if deletionProbability(jj)==0.25
%             legend('0%','5%','10%','15%','20%','25%');
%         end
        
        
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

%% Degradation results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Degradation results (RMSE): Frequency indexes (Lomb), random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i          ',deletionProbability(2:end)*100); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLomb{kk,jj}.lf; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltLomb{kk-length(subject),jj}.lf; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('P_{LF} [ms^2]')
% RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)); %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_lomb_' fillGaps '_lf']),'epsc'); end %#ok<*UNRCH>
fprintf('P_LF (norm)    '); fprintf('%.3f       ',nRMSE(2:end)); fprintf('\n')


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLomb{kk,jj}.hf; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltLomb{kk-length(subject),jj}.hf; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('P_{HF} [ms^2]')
% RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)); %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_lomb_' fillGaps '_hf']),'epsc'); end
fprintf('P_HF (norm)    '); fprintf('%.3f       ',nRMSE(2:end)); fprintf('\n')


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLomb{kk,jj}.lfn; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltLomb{kk-length(subject),jj}.lfn; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('P_{LFn}')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))); % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_lomb_' fillGaps '_lfn']),'epsc'); end
fprintf('P_LFn          '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLomb{kk,jj}.lfhf; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltLomb{kk-length(subject),jj}.lfhf; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('P_{LF}/P_{HF}')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))); % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_lomb_' fillGaps '_lfhf']),'epsc'); end
fprintf('P_LF/P_HF      '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')

%% Sympathovagal balance results

% aux1 = [];
% aux2 = [];
% significance = [];
% for jj=1:length(deletionProbability)
%     for kk=1:length(subject)
%         aux1(kk,jj) = resultsSupineLomb{kk,jj}.lf;
%         aux2(kk,jj) = resultsTiltLomb{kk,jj}.lf;
%     end
%     significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
% end
% fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
% ylabel('P_{LF} [ms^2]')
% if saveFigs, saveas(gcf,strcat(['images\ans_rd_lomb_' fillGaps '_lf']),'epsc'); end
% 
% 
% aux1 = [];
% aux2 = [];
% significance = [];
% for jj=1:length(deletionProbability)
%     for kk=1:length(subject)
%         aux1(kk,jj) = resultsSupineLomb{kk,jj}.hf;
%         aux2(kk,jj) = resultsTiltLomb{kk,jj}.hf;
%     end
%     significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
% end
% fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
% ylabel('P_{HF} [ms^2]')
% if saveFigs, saveas(gcf,strcat(['images\ans_rd_lomb_' fillGaps '_hf']),'epsc'); end
% 
% 
% aux1 = [];
% aux2 = [];
% significance = [];
% for jj=1:length(deletionProbability)
%     for kk=1:length(subject)
%         aux1(kk,jj) = resultsSupineLomb{kk,jj}.lfn;
%         aux2(kk,jj) = resultsTiltLomb{kk,jj}.lfn;
%     end
%     significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
% end
% fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
% ylabel('P_{LFn}')
% if saveFigs, saveas(gcf,strcat(['images\ans_rd_lomb_' fillGaps '_lfn']),'epsc'); end
% 
% 
% aux1 = [];
% aux2 = [];
% significance = [];
% for jj=1:length(deletionProbability)
%     for kk=1:length(subject)
%         aux1(kk,jj) = resultsSupineLomb{kk,jj}.lfhf; %#ok<*SAGROW>
%         aux2(kk,jj) = resultsTiltLomb{kk,jj}.lfhf;
%     end
%     significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
% end
% fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
% ylabel('P_{LF}/P_{HF}')
% if saveFigs, saveas(gcf,strcat(['images\ans_rd_lomb_' fillGaps '_lfhf']),'epsc'); end

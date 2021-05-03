%% AR results with random distributed errors

clear

% Signal duration
SIGNAL_DURATION = 120; % [s]

% subject = {'07M'};
subject = {'01M' '02V' '04V' '05V' '06M' '07M' '11M' '13V' '17V'}; % Respiración en banda de HF
% subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '11M'...
%     '12V' '13V' '15V' '16V' '17V'};
deletionProbability = [0 0.05 0.1 0.15 0.2 0.25];
fillGaps = 'iterative'; % 'incidences' 'iterative'
saveFigs = false;

resultsSupineAR = cell(length(subject),length(deletionProbability));
resultsTiltAR = cell(length(subject),length(deletionProbability));
for kk = 1:length(subject)
    for jj = 1:length(deletionProbability)
        load(strcat(pwd,'\Databases\TiltPreproc\',subject{kk},'.mat'));
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
                    [~,~,~,SupineECG.tn] = incidences(SupineECG.tk,1); SupineECG.ids = [];
                    [~,~,~,TiltECG.tn] = incidences(TiltECG.tk,1); TiltECG.ids = [];
            end
        else
            SupineECG.tn = SupineECG.tk; SupineECG.ids = [];
            TiltECG.tn = TiltECG.tk; TiltECG.ids = [];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % AR Model
        SupineAR = ARModelEstimation(diff(SupineECG.tn));
        TiltAR = ARModelEstimation(diff(TiltECG.tn));
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         if deletionProbability(jj)==0
%             figure('DefaultAxesFontSize',14); hold on;
%             plot(TiltAR.f,TiltAR.psd,'LineWidth',2); hold on
% %             title(subject{kk})
%         else
%             plot(TiltAR.f,TiltAR.psd,'LineWidth',1); hold on
%         end
%         xlabel('Frequency [Hz]')
%         ylabel('PSD_{AR} [ms^2/Hz]')
%         axis tight
%         xlim([0 0.4])
%         ylim([0 12000])
%         set(gcf,'position',[0,0,1000,300])
%         if deletionProbability(jj)==0.25
%             legend('0%','5%','10%','15%','20%','25%');
%         end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
        
        % Save results
        resultsSupineAR{kk,jj} = SupineAR;
        resultsTiltAR{kk,jj} = TiltAR;           
    end
end

%% Degradation results

aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineAR{kk,jj}.lf; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltAR{kk-length(subject),jj}.lf; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('P_{LF} [ms^2]')
% RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_ar_' fillGaps '_lf']),'epsc'); end %#ok<*UNRCH>


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineAR{kk,jj}.hf; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltAR{kk-length(subject),jj}.hf; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('P_{HF} [ms^2]')
% RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_ar_' fillGaps '_hf']),'epsc'); end


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineAR{kk,jj}.lfn; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltAR{kk-length(subject),jj}.lfn; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('P_{LFn}')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_ar_' fillGaps '_lfn']),'epsc'); end


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineAR{kk,jj}.lfhf; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltAR{kk-length(subject),jj}.lfhf; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('P_{LF}/P_{HF}')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_ar_' fillGaps '_lfhf']),'epsc'); end


%% Sympathovagal balance results

aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineAR{kk,jj}.lf;
        aux2(kk,jj) = resultsTiltAR{kk,jj}.lf;
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
end
fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
ylabel('P_{LF} [ms^2]')
if saveFigs, saveas(gcf,strcat(['images\ans_rd_ar_' fillGaps '_lf']),'epsc'); end


aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineAR{kk,jj}.hf;
        aux2(kk,jj) = resultsTiltAR{kk,jj}.hf;
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
end
fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
ylabel('P_{HF} [ms^2]')
if saveFigs, saveas(gcf,strcat(['images\ans_rd_ar_' fillGaps '_hf']),'epsc'); end


aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineAR{kk,jj}.lfn;
        aux2(kk,jj) = resultsTiltAR{kk,jj}.lfn;
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
end
fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
ylabel('P_{LFn}')
if saveFigs, saveas(gcf,strcat(['images\ans_rd_ar_' fillGaps '_lfn']),'epsc'); end


aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineAR{kk,jj}.lfhf; %#ok<*SAGROW>
        aux2(kk,jj) = resultsTiltAR{kk,jj}.lfhf;
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
end
fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
ylabel('P_{LF}/P_{HF}')
if saveFigs, saveas(gcf,strcat(['images\ans_rd_ar_' fillGaps '_lfhf']),'epsc'); end

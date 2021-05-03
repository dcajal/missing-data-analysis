%% Symbolic results with random distributed errors

clear

% Signal duration
SIGNAL_DURATION = 120; % [s]

subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '11M'...
    '12V' '13V' '15V' '16V' '17V'};
burstDuration = [0 10 20 30 40 50 60];
fillGaps = 'none'; % 'none' 'incidences' 'iterative'
saveFigs = false;

resultsSupineHRF = cell(length(subject),length(burstDuration));
resultsTiltHRF = cell(length(subject),length(burstDuration));
for kk = 1:length(subject)
    for jj = 1:length(burstDuration)
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

        if burstDuration(jj) > 0
            rng('default')
    
            SupineECG.tk(find(SupineECG.tk<=60-burstDuration(jj)/2,1,'Last'):...
                find(SupineECG.tk>=60+burstDuration(jj)/2,1)) = [];
            TiltECG.tk(find(TiltECG.tk<=60-burstDuration(jj)/2,1,'Last'):...
                find(TiltECG.tk>=60+burstDuration(jj)/2,1)) = [];

            switch fillGaps
                case 'iterative'
                    SupineECG.tn = gapcorrector(SupineECG.tk);
                    TiltECG.tn = gapcorrector(TiltECG.tk);
                case 'incidences'
                    [~,~,~,SupineECG.tn] = incidences(SupineECG.tk,1);
                    [~,~,~,TiltECG.tn] = incidences(TiltECG.tk,1);
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

        % Symbolic metrics
        SupineHRF = hrf(SupineECG.tn);
        TiltHRF = hrf(TiltECG.tn);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Save results
        resultsSupineHRF{kk,jj} = SupineHRF;
        resultsTiltHRF{kk,jj} = TiltHRF;   
    end
end

%% Degradation results

aux1 = [];
significance = [];
for jj=1:length(burstDuration)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineHRF{kk,jj}.pip; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltHRF{kk-length(subject),jj}.pip; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(burstDuration));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,burstDuration,significance)
ylabel('PIP [%]')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_bursts_hrf_' fillGaps '_pip']),'epsc'); end


aux1 = [];
significance = [];
for jj=1:length(burstDuration)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineHRF{kk,jj}.ials; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltHRF{kk-length(subject),jj}.ials; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(burstDuration));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,burstDuration,significance)
ylabel('IALS [intervals^{-1}]')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_bursts_hrf_' fillGaps '_ials']),'epsc'); end


aux1 = [];
significance = [];
for jj=1:length(burstDuration)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineHRF{kk,jj}.pss; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltHRF{kk-length(subject),jj}.pss; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(burstDuration));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,burstDuration,significance)
ylabel('PSS [%]')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_bursts_hrf_' fillGaps '_pss']),'epsc'); end


aux1 = [];
significance = [];
for jj=1:length(burstDuration)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineHRF{kk,jj}.pas; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltHRF{kk-length(subject),jj}.pas; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(burstDuration));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,burstDuration,significance)
ylabel('PAS [%]')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_bursts_hrf_' fillGaps '_pas']),'epsc'); end


%% Sympathovagal balance results

aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(burstDuration)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineHRF{kk,jj}.pip;
        aux2(kk,jj) = resultsTiltHRF{kk,jj}.pip;
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
end
fancyBoxplot(aux1,aux2,burstDuration,significance,true)
ylabel('PIP [%]')
if saveFigs, saveas(gcf,strcat(['images\ans_bursts_hrf_' fillGaps '_pip']),'epsc'); end


aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(burstDuration)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineHRF{kk,jj}.ials;
        aux2(kk,jj) = resultsTiltHRF{kk,jj}.ials;
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
end
fancyBoxplot(aux1,aux2,burstDuration,significance,true)
ylabel('IALS [intervals^{-1}]')
if saveFigs, saveas(gcf,strcat(['images\ans_bursts_hrf_' fillGaps '_ials']),'epsc'); end


aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(burstDuration)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineHRF{kk,jj}.pss;
        aux2(kk,jj) = resultsTiltHRF{kk,jj}.pss;
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
end
fancyBoxplot(aux1,aux2,burstDuration,significance,true)
ylabel('PSS [%]')
if saveFigs, saveas(gcf,strcat(['images\ans_bursts_hrf_' fillGaps '_pss']),'epsc'); end


aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(burstDuration)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineHRF{kk,jj}.pas;
        aux2(kk,jj) = resultsTiltHRF{kk,jj}.pas;
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
end
fancyBoxplot(aux1,aux2,burstDuration,significance,true)
ylabel('PAS [%]')
if saveFigs, saveas(gcf,strcat(['images\ans_bursts_hrf_' fillGaps '_pas']),'epsc'); end

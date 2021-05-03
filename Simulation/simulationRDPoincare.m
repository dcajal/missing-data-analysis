%% Poincare results with random distributed errors

clear

% Signal duration
SIGNAL_DURATION = 120; % [s]

subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '11M'...
    '12V' '13V' '15V' '16V' '17V'};
deletionProbability = [0 0.05 0.1 0.15 0.2 0.25];
fillGaps = 'none'; % 'incidences' 'iterative' 'none'
detectGaps = true;
saveFigs = false;

resultsSupineLPP = cell(length(subject),length(deletionProbability));
resultsTiltLPP = cell(length(subject),length(deletionProbability));
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
        
        % Lagged Poincare Plots
        SupineLPP = lpp(SupineECG.tn,1,detectGaps);
        TiltLPP = lpp(TiltECG.tn,1,detectGaps);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
        % Save results
        resultsSupineLPP{kk,jj} = SupineLPP;
        resultsTiltLPP{kk,jj} = TiltLPP;            
    end
end

%% Degradation results

aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLPP{kk,jj}.SD1; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltLPP{kk-length(subject),jj}.SD1; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('SD1 [ms]')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_lpp_' fillGaps '_sd1']),'epsc'); end %#ok<*UNRCH>


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLPP{kk,jj}.SD2; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltLPP{kk-length(subject),jj}.SD2; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('SD2 [ms]')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_lpp_' fillGaps '_sd2']),'epsc'); end


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLPP{kk,jj}.SD12; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltLPP{kk-length(subject),jj}.SD12; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('SD12')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_lpp_' fillGaps '_sd12']),'epsc'); end


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLPP{kk,jj}.S; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltLPP{kk-length(subject),jj}.S; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('S [ms^2]')
% RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_lpp_' fillGaps '_s']),'epsc'); end


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLPP{kk,jj}.Md; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltLPP{kk-length(subject),jj}.Md; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('Md [ms]')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_lpp_' fillGaps '_md']),'epsc'); end


aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLPP{kk,jj}.Sd; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltLPP{kk-length(subject),jj}.Sd; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('Sd [ms]')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))) % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_lpp_' fillGaps '_sd']),'epsc'); end


%% Sympathovagal balance results

aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLPP{kk,jj}.SD1;
        aux2(kk,jj) = resultsTiltLPP{kk,jj}.SD1;
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
end
fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
ylabel('SD1 [ms]')
if saveFigs, saveas(gcf,strcat(['images\ans_rd_lpp_' fillGaps '_sd1']),'epsc'); end


aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLPP{kk,jj}.SD2;
        aux2(kk,jj) = resultsTiltLPP{kk,jj}.SD2;
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
end
fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
ylabel('SD2 [ms]')
if saveFigs, saveas(gcf,strcat(['images\ans_rd_lpp_' fillGaps '_sd2']),'epsc'); end


aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLPP{kk,jj}.SD12;
        aux2(kk,jj) = resultsTiltLPP{kk,jj}.SD12;
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
end
fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
ylabel('SD12')
if saveFigs, saveas(gcf,strcat(['images\ans_rd_lpp_' fillGaps '_sd12']),'epsc'); end


aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLPP{kk,jj}.S;
        aux2(kk,jj) = resultsTiltLPP{kk,jj}.S;
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
end
fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
ylabel('S [ms^2]')
if saveFigs, saveas(gcf,strcat(['images\ans_rd_lpp_' fillGaps '_s']),'epsc'); end


aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLPP{kk,jj}.Md;
        aux2(kk,jj) = resultsTiltLPP{kk,jj}.Md;
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
end
fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
ylabel('Md [ms]')
if saveFigs, saveas(gcf,strcat(['images\ans_rd_lpp_' fillGaps '_md']),'epsc'); end


aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineLPP{kk,jj}.Sd;
        aux2(kk,jj) = resultsTiltLPP{kk,jj}.Sd;
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
end
fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
ylabel('Sd [ms]')
if saveFigs, saveas(gcf,strcat(['images\ans_rd_lpp_' fillGaps '_sd']),'epsc'); end


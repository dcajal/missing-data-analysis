%% Time-domain results with random distributed errors

clear
addpath('lib','database');
set(0,'defaultAxesFontName', 'Helvetica');
set(0,'defaultTextFontName', 'Helvetica');
set(0,'defaulttextInterpreter','latex')
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize', 20);
set(0,'defaultTextFontSize', 20);
set(0,'defaultFigureRenderer', 'painters')

% Signal duration
SIGNAL_DURATION = 120; % [s]

subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '11M'...
    '12V' '13V' '15V' '16V' '17V'};
deletionProbability = 0:0.05:0.25;
fillGaps = 'iterative'; % 'none' 'incidences' 'iterative'
detectGaps = false;
saveFigs = false;

resultsSupineTDP = cell(length(subject),length(deletionProbability));
resultsTiltTDP = cell(length(subject),length(deletionProbability));
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
                    SupineECG.tn = gapcorrector(SupineECG.tk,true);
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

        % Time domain parameters
        [SupineTDP] = tempind(SupineECG.tn,detectGaps);
        [TiltTDP] = tempind(TiltECG.tn,detectGaps);
        
        
        % Save results
        resultsSupineTDP{kk,jj} = SupineTDP;
        resultsTiltTDP{kk,jj} = TiltTDP;
             
    end
end

%% Degradation results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Degradation results (RMSE): Time indexes, random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i          ',deletionProbability(2:end)*100); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineTDP{kk,jj}.HRM; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltTDP{kk-length(subject),jj}.HRM; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('MHR [beats/min]','interpreter','tex')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))); % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_tdp_' fillGaps '_mhr']),'epsc'); end %#ok<*UNRCH>
fprintf('MHR              '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')

aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineTDP{kk,jj}.SDNN; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltTDP{kk-length(subject),jj}.SDNN; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('SDNN [ms]','interpreter','tex')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))); % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_tdp_' fillGaps '_sdnn']),'epsc'); end
fprintf('SDNN             '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')

aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineTDP{kk,jj}.SDSD; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltTDP{kk-length(subject),jj}.SDSD; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('SDSD [ms]','interpreter','tex')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))); % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_tdp_' fillGaps '_sdsd']),'epsc'); end
fprintf('SDSD             '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')

aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineTDP{kk,jj}.RMSSD; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltTDP{kk-length(subject),jj}.RMSSD; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('RMSSD [ms]','interpreter','tex')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))); % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_tdp_' fillGaps '_rmssd']),'epsc'); end
fprintf('RMSSD            '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')

aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:length(subject)
        aux1(kk,jj) = resultsSupineTDP{kk,jj}.pNN50; %#ok<*SAGROW>
    end
    for kk=length(subject)+1:2*length(subject)
        aux1(kk,jj) = resultsTiltTDP{kk-length(subject),jj}.pNN50; %#ok<*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj));
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
fancyBoxplot(aux2,aux1,100*deletionProbability,significance)
ylabel('pNN50 [%]','interpreter','tex')
RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject))); % RMSE
% nRMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*length(subject)))./nanmean(aux2(:,1)) %#ok<*NOPTS> % nRMSE
if saveFigs, saveas(gcf,strcat(['images\deg_rd_tdp_' fillGaps '_pnn50']),'epsc'); end
fprintf('pNN50            '); fprintf('%.3f       ',RMSE(2:end)); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')

%% Sympathovagal balance results
 
% fprintf('\n'); fprintf('\n')
% fprintf('---------------------------------------------------------------------------------------\n')
% fprintf('p-values (supine/tilt groups): Time indexes, random distributed errors, using %s method\n',fillGaps);
% disp('Measure          Deletion probability (%)');
% fprintf('                 '); fprintf('%i          ',deletionProbability(2:end)*100); fprintf('\n')
% fprintf('---------------------------------------------------------------------------------------\n')
% 
% aux1 = [];
% aux2 = [];
% significance = [];
% for jj=1:length(deletionProbability)
%     for kk=1:length(subject)
%         aux1(kk,jj) = resultsSupineTDP{kk,jj}.HRM;
%         aux2(kk,jj) = resultsTiltTDP{kk,jj}.HRM;
%     end
%     significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
% end
% fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
% ylabel('MHR [beats/min]','interpreter','tex')
% if saveFigs, saveas(gcf,strcat(['images\ans_rd_tdp_' fillGaps '_mhr']),'epsc'); end
% fprintf('MHR             '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')
% 
% aux1 = [];
% aux2 = [];
% significance = [];
% for jj=1:length(deletionProbability)
%     for kk=1:length(subject)
%         aux1(kk,jj) = resultsSupineTDP{kk,jj}.SDNN;
%         aux2(kk,jj) = resultsTiltTDP{kk,jj}.SDNN;
%     end
%     significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
% end
% fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
% ylabel('SDNN [ms]','interpreter','tex')
% if saveFigs, saveas(gcf,strcat(['images\ans_rd_tdp_' fillGaps '_sdnn']),'epsc'); end
% fprintf('SDNN            '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')
% 
% aux1 = [];
% aux2 = [];
% significance = [];
% for jj=1:length(deletionProbability)
%     for kk=1:length(subject)
%         aux1(kk,jj) = resultsSupineTDP{kk,jj}.SDSD;
%         aux2(kk,jj) = resultsTiltTDP{kk,jj}.SDSD;
%     end
%     significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
% end
% fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
% ylabel('SDSD [ms]','interpreter','tex')
% if saveFigs, saveas(gcf,strcat(['images\ans_rd_tdp_' fillGaps '_sdsd']),'epsc'); end
% fprintf('SDSD            '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')
% 
% aux1 = [];
% aux2 = [];
% significance = [];
% for jj=1:length(deletionProbability)
%     for kk=1:length(subject)
%         aux1(kk,jj) = resultsSupineTDP{kk,jj}.RMSSD;
%         aux2(kk,jj) = resultsTiltTDP{kk,jj}.RMSSD;
%     end
%     significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
% end
% fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
% ylabel('RMSSD [ms]','interpreter','tex')
% if saveFigs, saveas(gcf,strcat(['images\ans_rd_tdp_' fillGaps '_rmssd']),'epsc'); end
% fprintf('RMSSD           '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')
% 
% aux1 = [];
% aux2 = [];
% significance = [];
% for jj=1:length(deletionProbability)
%     for kk=1:length(subject)
%         aux1(kk,jj) = resultsSupineTDP{kk,jj}.pNN50;
%         aux2(kk,jj) = resultsTiltTDP{kk,jj}.pNN50;
%     end
%     significance(jj) = signrank(aux1(:,jj),aux2(:,jj));
% end
% fancyBoxplot(aux1,aux2,100*deletionProbability,significance,true)
% ylabel('pNN50 [\%]','interpreter','tex')
% if saveFigs, saveas(gcf,strcat(['images\ans_rd_tdp_' fillGaps '_pnn50']),'epsc'); end
% fprintf('pNN50           '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')


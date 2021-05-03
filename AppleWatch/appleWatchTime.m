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

debug = false;
segmentDurationSec = 2*60;
stepDurationSec = 1*60;

% Load all files in database directory
dirlist = dir('database');
files = cell([1 length(dirlist)-2]);
for kk = 3:length(dirlist)
    files{kk-2} = dirlist(kk).name;
end

saveIndex = 1;
for kk = 1:length(files)
    load(strcat(pwd,'\database\',files{kk}),'tRRAWrelax','tRRH7relax');
    tRRAWrelax = tRRAWrelax/1000;
    tRRH7relax = tRRH7relax/1000; 
    
    segmentBegin = 0:stepDurationSec:tRRH7relax(end)-segmentDurationSec;
    segmentEnd = segmentBegin+segmentDurationSec;
    
    if ~isempty(segmentBegin)
        for jj = 1:2%length(segmentBegin)
            fprintf('Analyzing %s (relax), Segment %d...',files{kk},jj);
            referenceSegment = tRRH7relax(tRRH7relax>segmentBegin(jj) & tRRH7relax<segmentEnd(jj));
            awSegment = tRRAWrelax(tRRAWrelax>segmentBegin(jj) & tRRAWrelax<segmentEnd(jj));
            
            % Reference
            reference{saveIndex} = tempind(referenceSegment,true); %#ok<SAGROW>
            
            % No preprocessing
            awNoPreproc{saveIndex} = tempind(awSegment,false); %#ok<SAGROW>
            
            % Detecting gaps
            awRemovingOutliers{saveIndex} = tempind(awSegment,true); %#ok<SAGROW>
            
            % Incidences method
            [~,~,~,qrsIncidences] = incidences(awSegment);
            awIncidences{saveIndex} = tempind(qrsIncidences,false); clear aux %#ok<SAGROW>
            
            % Iterative method
            qrsIterative = gapcorrector(awSegment);
            awIterative{saveIndex} = tempind(qrsIterative,true); %#ok<SAGROW>
%             
%             figure;
%             ax(1) = subplot(411); 
%             stem(referenceSegment(2:end), diff(referenceSegment)); ylabel('Holter')
%             ax(2) = subplot(412);
%             stem(awSegment(2:end), diff(awSegment)); ylabel('Armband')
%             ax(3) = subplot(413);
%             stem(qrsIncidences(2:end), diff(qrsIncidences)); ylabel('Incidences correction')
%             ax(4) = subplot(414);
%             stem(qrsIterative(2:end), diff(qrsIterative)); ylabel('Iterative correction')
%             linkaxes(ax,'x'); axis tight; set(gcf, 'Position', get(0, 'Screensize'));
%             pause
                        
            saveIndex = saveIndex+1;
            fprintf('Done\n');
        end
    end
end

for kk = 1:length(files)
    load(strcat(pwd,'\database\',files{kk}),'tRRAWstress','tRRH7stress');  
    tRRAWstress = tRRAWstress/1000;
    tRRH7stress = tRRH7stress/1000;  
    
    segmentBegin = 0:stepDurationSec:tRRH7stress(end)-segmentDurationSec;
    segmentEnd = segmentBegin+segmentDurationSec;
    
    if ~isempty(segmentBegin)
        for jj = 1:2%length(segmentBegin)
            fprintf('Analyzing %s (stress), Segment %d...',files{kk},jj);
            referenceSegment = tRRH7stress(tRRH7stress>segmentBegin(jj) & tRRH7stress<segmentEnd(jj));
            awSegment = tRRAWstress(tRRAWstress>segmentBegin(jj) & tRRAWstress<segmentEnd(jj));
            
            % Reference
            reference{saveIndex} = tempind(referenceSegment,true);
            
            % No preprocessing
            awNoPreproc{saveIndex} = tempind(awSegment,false);
            
            % Detecting gaps
            awRemovingOutliers{saveIndex} = tempind(awSegment,true);
            
            % Incidences method
            [~,~,~,qrsIncidences] = incidences(awSegment);
            awIncidences{saveIndex} = tempind(qrsIncidences,false); clear aux
            
            % Iterative method
            qrsIterative = gapcorrector(awSegment);
            awIterative{saveIndex} = tempind(qrsIterative,true);
%                 
%             figure;
%             ax(1) = subplot(411);
%             stem(referenceSegment(2:end), diff(referenceSegment)); ylabel('Holter')
%             ax(2) = subplot(412);
%             stem(awSegment(2:end), diff(awSegment)); ylabel('Armband')
%             ax(3) = subplot(413);
%             stem(qrsIncidences(2:end), diff(qrsIncidences)); ylabel('Incidences correction')
%             ax(4) = subplot(414);
%             stem(qrsIterative(2:end), diff(qrsIterative)); ylabel('Iterative correction')
%             linkaxes(ax,'x'); axis tight;
%             pause

            saveIndex = saveIndex+1;
            fprintf('Done\n');
        end
    end
end

%%
errorThreshold = 0.001:0.001:0.1;
% errorThreshold = 0.01;

%% MHR

noPreprocCorrects = zeros(length(reference),numel(errorThreshold));
removingOutliersCorrects = zeros(length(reference),numel(errorThreshold));
incidencesCorrects = zeros(length(reference),numel(errorThreshold));
iterativeCorrects = zeros(length(reference),numel(errorThreshold));
for kk=1:length(reference)
    noPreprocCorrects(kk,:) = abs((awNoPreproc{kk}.HRM - reference{kk}.HRM)/reference{kk}.HRM)<errorThreshold;
    removingOutliersCorrects(kk,:) = abs((awRemovingOutliers{kk}.HRM - reference{kk}.HRM)/reference{kk}.HRM)<errorThreshold;
    incidencesCorrects(kk,:) = abs((awIncidences{kk}.HRM - reference{kk}.HRM)/reference{kk}.HRM)<errorThreshold;
    iterativeCorrects(kk,:) = abs((awIterative{kk}.HRM - reference{kk}.HRM)/reference{kk}.HRM)<errorThreshold;
end

% fprintf('\n');
% fprintf('----------------------------------------\n');
% fprintf('---------------- MHR -------------------\n');
% fprintf('----------------------------------------\n');
% fprintf('No preprocessing coverage: %i out of %i (%.2f%%)\n',sum(noPreprocCorrects),numel(noPreprocCorrects),sum(noPreprocCorrects)/numel(noPreprocCorrects)*100);
% fprintf('Removing outliers coverage: %i out of %i (%.2f%%)\n',sum(removingOutliersCorrects),numel(removingOutliersCorrects),sum(removingOutliersCorrects)/numel(removingOutliersCorrects)*100);
% fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
% fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
% fprintf('----------------------------------------\n');
% fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
% fprintf('Iterative coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(iterativeCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
% fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
% fprintf('Incidences coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(incidencesCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
% fprintf('----------------------------------------\n');
% 
% figure; hold on
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 4 4],[1 1 1]);
% p(1) = bar(1:length(reference),noPreprocCorrects*4,'b');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 3 3],[1 1 1]);
% p(2) = bar(1:length(reference),removingOutliersCorrects*3,'g');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 2 2],[1 1 1]);
% p(3) = bar(1:length(reference),incidencesCorrects*2,'k');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 1 1],[1 1 1]);
% p(4) = bar(1:length(reference),iterativeCorrects,'r');
% legend(p,'No Prep','Removing Outliers','Incidences','Iterative','Location','bestoutside')
% title('MHR'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

figure; hold on
plot(errorThreshold*100,sum(noPreprocCorrects,1),'-');
plot(errorThreshold*100,sum(removingOutliersCorrects,1),'--');
plot(errorThreshold*100,sum(incidencesCorrects,1),'-*');
plot(errorThreshold*100,sum(iterativeCorrects,1),'-^');
legend('No Prep','Removing Outliers','Incidences','Iterative','Location','best')
title('MHR'); xlabel('Permited error (\%)'); ylabel('Number of correct cases')



%% SDNN

noPreprocCorrects = zeros(length(reference),numel(errorThreshold));
removingOutliersCorrects = zeros(length(reference),numel(errorThreshold));
incidencesCorrects = zeros(length(reference),numel(errorThreshold));
iterativeCorrects = zeros(length(reference),numel(errorThreshold));
for kk=1:length(reference)
    noPreprocCorrects(kk,:) = abs((awNoPreproc{kk}.SDNN - reference{kk}.SDNN)/reference{kk}.SDNN)<errorThreshold;
    removingOutliersCorrects(kk,:) = abs((awRemovingOutliers{kk}.SDNN - reference{kk}.SDNN)/reference{kk}.SDNN)<errorThreshold;
    incidencesCorrects(kk,:) = abs((awIncidences{kk}.SDNN - reference{kk}.SDNN)/reference{kk}.SDNN)<errorThreshold;
    iterativeCorrects(kk,:) = abs((awIterative{kk}.SDNN - reference{kk}.SDNN)/reference{kk}.SDNN)<errorThreshold;
end

% fprintf('\n');
% fprintf('----------------------------------------\n');
% fprintf('--------------- SDNN -------------------\n');
% fprintf('No preprocessing coverage: %i out of %i (%.2f%%)\n',sum(noPreprocCorrects),numel(noPreprocCorrects),sum(noPreprocCorrects)/numel(noPreprocCorrects)*100);
% fprintf('Removing outliers coverage: %i out of %i (%.2f%%)\n',sum(removingOutliersCorrects),numel(removingOutliersCorrects),sum(removingOutliersCorrects)/numel(removingOutliersCorrects)*100);
% fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
% fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
% fprintf('----------------------------------------\n');
% fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
% fprintf('Iterative coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(iterativeCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
% fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
% fprintf('Incidences coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(incidencesCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
% fprintf('----------------------------------------\n');

% figure; hold on
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 4 4],[1 1 1]);
% p(1) = bar(1:length(reference),noPreprocCorrects*4,'b');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 3 3],[1 1 1]);
% p(2) = bar(1:length(reference),removingOutliersCorrects*3,'g');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 2 2],[1 1 1]);
% p(3) = bar(1:length(reference),incidencesCorrects*2,'k');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 1 1],[1 1 1]);
% p(4) = bar(1:length(reference),iterativeCorrects,'r');
% legend(p,'No Prep','Removing Outliers','Incidences','Iterative','Location','bestoutside')
% title('SDNN'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

figure; hold on
plot(errorThreshold*100,sum(noPreprocCorrects,1),'-');
plot(errorThreshold*100,sum(removingOutliersCorrects,1),'--');
plot(errorThreshold*100,sum(incidencesCorrects,1),'-*');
plot(errorThreshold*100,sum(iterativeCorrects,1),'-^');
legend('No Prep','Removing Outliers','Incidences','Iterative','Location','best')
title('SDNN'); xlabel('Permited error (\%)'); ylabel('Number of correct cases')

%% RMSSD

noPreprocCorrects = zeros(length(reference),numel(errorThreshold));
removingOutliersCorrects = zeros(length(reference),numel(errorThreshold));
incidencesCorrects = zeros(length(reference),numel(errorThreshold));
iterativeCorrects = zeros(length(reference),numel(errorThreshold));
for kk=1:length(reference)
    noPreprocCorrects(kk,:) = abs((awNoPreproc{kk}.RMSSD - reference{kk}.RMSSD)/reference{kk}.RMSSD)<errorThreshold;
    removingOutliersCorrects(kk,:) = abs((awRemovingOutliers{kk}.RMSSD - reference{kk}.RMSSD)/reference{kk}.RMSSD)<errorThreshold;
    incidencesCorrects(kk,:) = abs((awIncidences{kk}.RMSSD - reference{kk}.RMSSD)/reference{kk}.RMSSD)<errorThreshold;
    iterativeCorrects(kk,:) = abs((awIterative{kk}.RMSSD - reference{kk}.RMSSD)/reference{kk}.RMSSD)<errorThreshold;
end

% fprintf('\n');
% fprintf('----------------------------------------\n');
% fprintf('--------------- RMSSD ------------------\n');
% fprintf('----------------------------------------\n');
% fprintf('No preprocessing coverage: %i out of %i (%.2f%%)\n',sum(noPreprocCorrects),numel(noPreprocCorrects),sum(noPreprocCorrects)/numel(noPreprocCorrects)*100);
% fprintf('Removing outliers coverage: %i out of %i (%.2f%%)\n',sum(removingOutliersCorrects),numel(removingOutliersCorrects),sum(removingOutliersCorrects)/numel(removingOutliersCorrects)*100);
% fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
% fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
% fprintf('----------------------------------------\n');
% fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
% fprintf('Iterative coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(iterativeCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
% fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
% fprintf('Incidences coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(incidencesCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
% fprintf('----------------------------------------\n');

% figure; hold on
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 4 4],[1 1 1]);
% p(1) = bar(1:length(reference),noPreprocCorrects*4,'b');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 3 3],[1 1 1]);
% p(2) = bar(1:length(reference),removingOutliersCorrects*3,'g');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 2 2],[1 1 1]);
% p(3) = bar(1:length(reference),incidencesCorrects*2,'k');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 1 1],[1 1 1]);
% p(4) = bar(1:length(reference),iterativeCorrects,'r');
% legend(p,'No Prep','Removing Outliers','Incidences','Iterative','Location','bestoutside')
% title('RMSSD'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

figure; hold on
plot(errorThreshold*100,sum(noPreprocCorrects,1),'-');
plot(errorThreshold*100,sum(removingOutliersCorrects,1),'--');
plot(errorThreshold*100,sum(incidencesCorrects,1),'-*');
plot(errorThreshold*100,sum(iterativeCorrects,1),'-^');
legend('No Prep','Removing Outliers','Incidences','Iterative','Location','best')
title('RMSSD'); xlabel('Permited error (\%)'); ylabel('Number of correct cases')

%% pNN50

noPreprocCorrects = zeros(length(reference),numel(errorThreshold));
removingOutliersCorrects = zeros(length(reference),numel(errorThreshold));
incidencesCorrects = zeros(length(reference),numel(errorThreshold));
iterativeCorrects = zeros(length(reference),numel(errorThreshold));
for kk=1:length(reference)
    noPreprocCorrects(kk,:) = abs((awNoPreproc{kk}.pNN50 - reference{kk}.pNN50)/reference{kk}.pNN50)<errorThreshold;
    removingOutliersCorrects(kk,:) = abs((awRemovingOutliers{kk}.pNN50 - reference{kk}.pNN50)/reference{kk}.pNN50)<errorThreshold;
    incidencesCorrects(kk,:) = abs((awIncidences{kk}.pNN50 - reference{kk}.pNN50)/reference{kk}.pNN50)<errorThreshold;
    iterativeCorrects(kk,:) = abs((awIterative{kk}.pNN50 - reference{kk}.pNN50)/reference{kk}.pNN50)<errorThreshold;
end

% fprintf('\n');
% fprintf('----------------------------------------\n');
% fprintf('--------------- pNN50 ------------------\n');
% fprintf('----------------------------------------\n');
% fprintf('No preprocessing coverage: %i out of %i (%.2f%%)\n',sum(noPreprocCorrects),numel(noPreprocCorrects),sum(noPreprocCorrects)/numel(noPreprocCorrects)*100);
% fprintf('Removing outliers coverage: %i out of %i (%.2f%%)\n',sum(removingOutliersCorrects),numel(removingOutliersCorrects),sum(removingOutliersCorrects)/numel(removingOutliersCorrects)*100);
% fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
% fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
% fprintf('----------------------------------------\n');
% fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
% fprintf('Iterative coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(iterativeCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
% fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
% fprintf('Incidences coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(incidencesCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
% fprintf('----------------------------------------\n');

% figure; hold on
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 4 4],[1 1 1]);
% p(1) = bar(1:length(reference),noPreprocCorrects*4,'b');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 3 3],[1 1 1]);
% p(2) = bar(1:length(reference),removingOutliersCorrects*3,'g');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 2 2],[1 1 1]);
% p(3) = bar(1:length(reference),incidencesCorrects*2,'k');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 1 1],[1 1 1]);
% p(4) = bar(1:length(reference),iterativeCorrects,'r');
% legend(p,'No Prep','Removing Outliers','Incidences','Iterative','Location','bestoutside')
% title('pNN50'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

figure; hold on
plot(100*errorThreshold,sum(noPreprocCorrects,1),'-');
plot(100*errorThreshold,sum(removingOutliersCorrects,1),'--');
plot(100*errorThreshold,sum(incidencesCorrects,1),'-*');
plot(100*errorThreshold,sum(iterativeCorrects,1),'-^');
legend('No Prep','Removing Outliers','Incidences','Iterative','Location','best')
title('pNN50'); xlabel('Permited error (\%)'); ylabel('Number of correct cases')



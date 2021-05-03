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
segmentDurationSec = 3*60;
stepDurationSec = 1*60;
errorThreshold = 0.1;

% Load all files in database directory
dirlist = dir('database');
files = cell([1 length(dirlist)-2]);
for kk = 3:length(dirlist)
    files{kk-2} = dirlist(kk).name;
end

saveIndex = 1;
for kk = 1:length(files)
    load(strcat(pwd,'\database\',files{kk}),'QRS1','QRS2','delay','art');
    QRS2 = QRS2+delay; clear delay
    
    segmentBegin = 0:stepDurationSec:QRS2(end)-segmentDurationSec;
    segmentEnd = segmentBegin+segmentDurationSec;
    
    % Count number of artifacts per segment
    artifactCounter = zeros(size(segmentBegin));
    for jj = 1:length(art.flag)
        if art.flag(jj)
            artifactCounter = artifactCounter+(segmentEnd>art.t(jj)-5 & segmentBegin<art.t(jj)+5);
        end
    end
    
    % Choose clean segments
    segmentBegin(artifactCounter>0) = [];
    segmentEnd(artifactCounter>0) = [];
    
    if ~isempty(segmentBegin)
        for jj = 1:length(segmentBegin)
            holterSegment = QRS2(QRS2>segmentBegin(jj) & QRS2<segmentEnd(jj));
            armbandSegment = QRS1(QRS1>segmentBegin(jj) & QRS1<segmentEnd(jj));
            
            % Reference
            holter{saveIndex} = tempind(holterSegment,true); %#ok<SAGROW>
            
            % No preprocessing
            armbandNoPreproc{saveIndex} = tempind(armbandSegment,false); %#ok<SAGROW>
            
            % Detecting gaps
            armbandRemovingOutliers{saveIndex} = tempind(armbandSegment,true); %#ok<SAGROW>
            
            % Incidences method
            [~,~,~,qrsIncidences] = incidences(armbandSegment);
            armbandIncidences{saveIndex} = tempind(qrsIncidences,false); clear aux %#ok<SAGROW>
            
            % Iterative method
            qrsIterative = gapcorrector(armbandSegment);
            armbandIterative{saveIndex} = tempind(qrsIterative,true); %#ok<SAGROW>
            
            if debug
                figure;
                ax(1) = subplot(411);
                stem(holterSegment(2:end), diff(holterSegment)); ylabel('Holter')
                ax(2) = subplot(412);
                stem(armbandSegment(2:end), diff(armbandSegment)); ylabel('Armband')
                ax(3) = subplot(413);
                stem(qrsIncidences(2:end), diff(qrsIncidences)); ylabel('Incidences correction')
                ax(4) = subplot(414);
                stem(qrsIterative(2:end), diff(qrsIterative)); ylabel('Iterative correction')
                linkaxes(ax,'x'); axis tight;
            end
            
            saveIndex = saveIndex+1;
        end
    end
end

%% MHR

noPreprocCorrects = zeros(length(holter),1);
removingOutliersCorrects = zeros(length(holter),1);
incidencesCorrects = zeros(length(holter),1);
iterativeCorrects = zeros(length(holter),1);
for kk=1:length(holter)
    noPreprocCorrects(kk) = abs((armbandNoPreproc{kk}.HRM - holter{kk}.HRM)/holter{kk}.HRM)<errorThreshold;
    removingOutliersCorrects(kk) = abs((armbandRemovingOutliers{kk}.HRM - holter{kk}.HRM)/holter{kk}.HRM)<errorThreshold;
    incidencesCorrects(kk) = abs((armbandIncidences{kk}.HRM - holter{kk}.HRM)/holter{kk}.HRM)<errorThreshold;
    iterativeCorrects(kk) = abs((armbandIterative{kk}.HRM - holter{kk}.HRM)/holter{kk}.HRM)<errorThreshold;
end
figure; hold on
p(1) = bar(1:length(holter),noPreprocCorrects*4,'b');
patch([1 length(holter) length(holter) 1],[0 0 3 3],[1 1 1]);
p(2) = bar(1:length(holter),removingOutliersCorrects*3,'g');
patch([1 length(holter) length(holter) 1],[0 0 2 2],[1 1 1]);
p(3) = bar(1:length(holter),incidencesCorrects*2,'k');
patch([1 length(holter) length(holter) 1],[0 0 1 1],[1 1 1]);
p(4) = bar(1:length(holter),iterativeCorrects,'r');
legend(p,'No Prep','Removing Outliers','Incidences','Iterative','Location','bestoutside')
title('MHR'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

fprintf('\n');
fprintf('----------------------------------------\n');
fprintf('---------------- MHR -------------------\n');
fprintf('----------------------------------------\n');
fprintf('No preprocessing coverage: %i out of %i (%.2f%%)\n',sum(noPreprocCorrects),numel(noPreprocCorrects),sum(noPreprocCorrects)/numel(noPreprocCorrects)*100);
fprintf('Removing outliers coverage: %i out of %i (%.2f%%)\n',sum(removingOutliersCorrects),numel(removingOutliersCorrects),sum(removingOutliersCorrects)/numel(removingOutliersCorrects)*100);
fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
fprintf('----------------------------------------\n');
fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
fprintf('Iterative coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(iterativeCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
fprintf('Incidences coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(incidencesCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
fprintf('----------------------------------------\n');


%% SDNN

noPreprocCorrects = zeros(length(holter),1);
removingOutliersCorrects = zeros(length(holter),1);
incidencesCorrects = zeros(length(holter),1);
iterativeCorrects = zeros(length(holter),1);
for kk=1:length(holter)
    noPreprocCorrects(kk) = abs((armbandNoPreproc{kk}.SDNN - holter{kk}.SDNN)/holter{kk}.SDNN)<errorThreshold;
    removingOutliersCorrects(kk) = abs((armbandRemovingOutliers{kk}.SDNN - holter{kk}.SDNN)/holter{kk}.SDNN)<errorThreshold;
    incidencesCorrects(kk) = abs((armbandIncidences{kk}.SDNN - holter{kk}.SDNN)/holter{kk}.SDNN)<errorThreshold;
    iterativeCorrects(kk) = abs((armbandIterative{kk}.SDNN - holter{kk}.SDNN)/holter{kk}.SDNN)<errorThreshold;
end
figure; hold on
p(1) = bar(1:length(holter),noPreprocCorrects*4,'b');
patch([1 length(holter) length(holter) 1],[0 0 3 3],[1 1 1]);
p(2) = bar(1:length(holter),removingOutliersCorrects*3,'g');
patch([1 length(holter) length(holter) 1],[0 0 2 2],[1 1 1]);
p(3) = bar(1:length(holter),incidencesCorrects*2,'k');
patch([1 length(holter) length(holter) 1],[0 0 1 1],[1 1 1]);
p(4) = bar(1:length(holter),iterativeCorrects,'r');
legend(p,'No Prep','Removing Outliers','Incidences','Iterative','Location','bestoutside')
title('SDNN'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

fprintf('\n');
fprintf('----------------------------------------\n');
fprintf('--------------- SDNN -------------------\n');
fprintf('No preprocessing coverage: %i out of %i (%.2f%%)\n',sum(noPreprocCorrects),numel(noPreprocCorrects),sum(noPreprocCorrects)/numel(noPreprocCorrects)*100);
fprintf('Removing outliers coverage: %i out of %i (%.2f%%)\n',sum(removingOutliersCorrects),numel(removingOutliersCorrects),sum(removingOutliersCorrects)/numel(removingOutliersCorrects)*100);
fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
fprintf('----------------------------------------\n');
fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
fprintf('Iterative coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(iterativeCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
fprintf('Incidences coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(incidencesCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
fprintf('----------------------------------------\n');

%% RMSSD

noPreprocCorrects = zeros(length(holter),1);
removingOutliersCorrects = zeros(length(holter),1);
incidencesCorrects = zeros(length(holter),1);
iterativeCorrects = zeros(length(holter),1);
for kk=1:length(holter)
    noPreprocCorrects(kk) = abs((armbandNoPreproc{kk}.RMSSD - holter{kk}.RMSSD)/holter{kk}.RMSSD)<errorThreshold;
    removingOutliersCorrects(kk) = abs((armbandRemovingOutliers{kk}.RMSSD - holter{kk}.RMSSD)/holter{kk}.RMSSD)<errorThreshold;
    incidencesCorrects(kk) = abs((armbandIncidences{kk}.RMSSD - holter{kk}.RMSSD)/holter{kk}.RMSSD)<errorThreshold;
    iterativeCorrects(kk) = abs((armbandIterative{kk}.RMSSD - holter{kk}.RMSSD)/holter{kk}.RMSSD)<errorThreshold;
end
figure; hold on
p(1) = bar(1:length(holter),noPreprocCorrects*4,'b');
patch([1 length(holter) length(holter) 1],[0 0 3 3],[1 1 1]);
p(2) = bar(1:length(holter),removingOutliersCorrects*3,'g');
patch([1 length(holter) length(holter) 1],[0 0 2 2],[1 1 1]);
p(3) = bar(1:length(holter),incidencesCorrects*2,'k');
patch([1 length(holter) length(holter) 1],[0 0 1 1],[1 1 1]);
p(4) = bar(1:length(holter),iterativeCorrects,'r');
legend(p,'No Prep','Removing Outliers','Incidences','Iterative','Location','bestoutside')
title('RMSSD'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

fprintf('\n');
fprintf('----------------------------------------\n');
fprintf('--------------- RMSSD ------------------\n');
fprintf('----------------------------------------\n');
fprintf('No preprocessing coverage: %i out of %i (%.2f%%)\n',sum(noPreprocCorrects),numel(noPreprocCorrects),sum(noPreprocCorrects)/numel(noPreprocCorrects)*100);
fprintf('Removing outliers coverage: %i out of %i (%.2f%%)\n',sum(removingOutliersCorrects),numel(removingOutliersCorrects),sum(removingOutliersCorrects)/numel(removingOutliersCorrects)*100);
fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
fprintf('----------------------------------------\n');
fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
fprintf('Iterative coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(iterativeCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
fprintf('Incidences coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(incidencesCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
fprintf('----------------------------------------\n');

%% pNN50

noPreprocCorrects = zeros(length(holter),1);
removingOutliersCorrects = zeros(length(holter),1);
incidencesCorrects = zeros(length(holter),1);
iterativeCorrects = zeros(length(holter),1);
for kk=1:length(holter)
    noPreprocCorrects(kk) = abs((armbandNoPreproc{kk}.pNN50 - holter{kk}.pNN50)/holter{kk}.pNN50)<errorThreshold;
    removingOutliersCorrects(kk) = abs((armbandRemovingOutliers{kk}.pNN50 - holter{kk}.pNN50)/holter{kk}.pNN50)<errorThreshold;
    incidencesCorrects(kk) = abs((armbandIncidences{kk}.pNN50 - holter{kk}.pNN50)/holter{kk}.pNN50)<errorThreshold;
    iterativeCorrects(kk) = abs((armbandIterative{kk}.pNN50 - holter{kk}.pNN50)/holter{kk}.pNN50)<errorThreshold;
end
figure; hold on
p(1) = bar(1:length(holter),noPreprocCorrects*4,'b');
patch([1 length(holter) length(holter) 1],[0 0 3 3],[1 1 1]);
p(2) = bar(1:length(holter),removingOutliersCorrects*3,'g');
patch([1 length(holter) length(holter) 1],[0 0 2 2],[1 1 1]);
p(3) = bar(1:length(holter),incidencesCorrects*2,'k');
patch([1 length(holter) length(holter) 1],[0 0 1 1],[1 1 1]);
p(4) = bar(1:length(holter),iterativeCorrects,'r');
legend(p,'No Prep','Removing Outliers','Incidences','Iterative','Location','bestoutside')
title('pNN50'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

fprintf('\n');
fprintf('----------------------------------------\n');
fprintf('--------------- pNN50 ------------------\n');
fprintf('----------------------------------------\n');
fprintf('No preprocessing coverage: %i out of %i (%.2f%%)\n',sum(noPreprocCorrects),numel(noPreprocCorrects),sum(noPreprocCorrects)/numel(noPreprocCorrects)*100);
fprintf('Removing outliers coverage: %i out of %i (%.2f%%)\n',sum(removingOutliersCorrects),numel(removingOutliersCorrects),sum(removingOutliersCorrects)/numel(removingOutliersCorrects)*100);
fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
fprintf('----------------------------------------\n');
fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
fprintf('Iterative coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(iterativeCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
fprintf('Incidences coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~removingOutliersCorrects)),sum(~removingOutliersCorrects),sum(incidencesCorrects(~removingOutliersCorrects))/sum(~removingOutliersCorrects)*100);
fprintf('----------------------------------------\n');


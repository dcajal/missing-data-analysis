%% Time-domain results with random distributed errors

clear
addpath('lib','AppleWatch/database');
figurePresets

segmentDurationSec = 2*60;
stepDurationSec = 1*60;

% Load all files in database directory
dirlist = dir('AppleWatch/database');
files = cell([1 length(dirlist)-2]);
for kk = 3:length(dirlist)
    files{kk-2} = dirlist(kk).name;
end

saveIndex = 1;
for kk = 1:length(files)
    load(strcat(pwd,'/AppleWatch/database/',files{kk}),'tRRAWrelax','tRRH7relax','tRRAWstress','tRRH7stress');
    tRRAWrelax = tRRAWrelax/1000;
    tRRH7relax = tRRH7relax/1000; 
    tRRAWstress = tRRAWstress/1000;
    tRRH7stress = tRRH7stress/1000;  
    
    segmentBegin = 0:stepDurationSec:min(tRRH7relax(end),tRRH7stress(end))-segmentDurationSec;
    segmentEnd = segmentBegin+segmentDurationSec;
    
    if ~isempty(segmentBegin)
        for jj = 1:2%length(segmentBegin)
            fprintf('Analyzing %s (relax), Segment %d...',files{kk},jj);
            referenceSegment = tRRH7relax(tRRH7relax>segmentBegin(jj) & tRRH7relax<segmentEnd(jj));
            awSegment = tRRAWrelax(tRRAWrelax>segmentBegin(jj) & tRRAWrelax<segmentEnd(jj));
            
            % Reference
            reference{saveIndex} = tempind(referenceSegment); %#ok<*SAGROW>
            
            % No preprocessing
            awNoPreproc{saveIndex} = tempind(awSegment);
            
            % Detecting gaps
            awRemovingOutliers{saveIndex} = tempind(awSegment,true);
            
            % Incidences method
            [~,~,~,qrsIncidences] = incidences(awSegment);
            awIncidences{saveIndex} = tempind(qrsIncidences); clear aux
            
            % Iterative method (linear)
            qrsIterative = gapcorrector(awSegment);
            awIterative{saveIndex} = tempind(qrsIterative);
            
            % Iterative method (non linear)
            qrsIterativeNL = gapcorrectorNonLinear(awSegment);
            awIterativeNL{saveIndex} = tempind(qrsIterativeNL);
             
                        
            saveIndex = saveIndex+1;
            fprintf('Done\n');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf('Analyzing %s (stress), Segment %d...',files{kk},jj);
            referenceSegment = tRRH7stress(tRRH7stress>segmentBegin(jj) & tRRH7stress<segmentEnd(jj));
            awSegment = tRRAWstress(tRRAWstress>segmentBegin(jj) & tRRAWstress<segmentEnd(jj));
            
            % Reference
            reference{saveIndex} = tempind(referenceSegment); %#ok<*SAGROW>
            
            % No preprocessing
            awNoPreproc{saveIndex} = tempind(awSegment);
            
            % Detecting gaps
            awRemovingOutliers{saveIndex} = tempind(awSegment,true);
            
            % Incidences method
            [~,~,~,qrsIncidences] = incidences(awSegment);
            awIncidences{saveIndex} = tempind(qrsIncidences); clear aux
            
            % Iterative method (linear)
            qrsIterative = gapcorrector(awSegment);
            awIterative{saveIndex} = tempind(qrsIterative);
            
            % Iterative method (non linear)
            qrsIterativeNL = gapcorrectorNonLinear(awSegment);
            awIterativeNL{saveIndex} = tempind(qrsIterativeNL);
                 

            saveIndex = saveIndex+1;
            fprintf('Done\n');
        end
    end
end

results = [reference; awNoPreproc; awRemovingOutliers; awIncidences; awIterative; awIterativeNL];


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


%% Results

errorThreshold = 0.001:0.001:0.25;
% errorThreshold = 0.01;

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Absolute error: Time indexes, AW\n');
disp('Measure          Method');
fprintf('                 No preproc            Removing outliers            Incidences            Iterative            Iterative NL\n');
fprintf('---------------------------------------------------------------------------------------\n')

error = computeTimeCoverage(results, errorThreshold, 'MHR');
fprintf('MHR              '); fprintf('%.2f (%.2f-%.2f)        ',error); fprintf('\n')

error = computeTimeCoverage(results, errorThreshold, 'SDNN');
fprintf('SDNN             '); fprintf('%.2f (%.2f-%.2f)        ',error); fprintf('\n')

error = computeTimeCoverage(results, errorThreshold, 'SDSD');
fprintf('SDSD             '); fprintf('%.2f (%.2f-%.2f)        ',error); fprintf('\n')

% error = computeTimeCoverage(results, errorThreshold, 'RMSSD');
% fprintf('RMSSD            '); fprintf('%.2f (%.2f-%.2f)        ',error); fprintf('\n')

error = computeTimeCoverage(results, errorThreshold, 'pNN50');
fprintf('pNN50            '); fprintf('%.2f (%.2f-%.2f)        ',error); fprintf('\n')



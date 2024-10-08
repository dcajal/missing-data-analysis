%% Time-domain results with random distributed errors

clear
addpath('lib','AppleWatch/database');
figurePresets

segmentDurationSec = 120;
stepDurationSec = 50;

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
        for jj = 1:3%length(segmentBegin)
            fprintf('Analyzing %s (relax), Segment %d...',files{kk},jj);
            referenceSegment = tRRH7relax(tRRH7relax>segmentBegin(jj) & tRRH7relax<segmentEnd(jj));
            awSegment = tRRAWrelax(tRRAWrelax>segmentBegin(jj) & tRRAWrelax<segmentEnd(jj));
            
            % Reference
            reference{saveIndex} = tempind(referenceSegment); %#ok<*SAGROW>
                       
            % OR
            OR{saveIndex} = tempind(awSegment,true);
            
            % L
            L{saveIndex} = tempind(gapcorrector(awSegment));
            
            % NL
            NL{saveIndex} = tempind(gapcorrectorNonLinear(awSegment));
             
                        
            saveIndex = saveIndex+1;
            fprintf('Done\n');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf('Analyzing %s (stress), Segment %d...',files{kk},jj);
            referenceSegment = tRRH7stress(tRRH7stress>segmentBegin(jj) & tRRH7stress<segmentEnd(jj));
            awSegment = tRRAWstress(tRRAWstress>segmentBegin(jj) & tRRAWstress<segmentEnd(jj));
            
            % Reference
            reference{saveIndex} = tempind(referenceSegment); %#ok<*SAGROW>
            
            % OR
            OR{saveIndex} = tempind(awSegment,true);
            
            % L
            L{saveIndex} = tempind(gapcorrector(awSegment));
            
            % NL
            NL{saveIndex} = tempind(gapcorrectorNonLinear(awSegment));
                 

            saveIndex = saveIndex+1;
            fprintf('Done\n');
        end
    end
end

results = [reference; OR; L; NL];


%% Reference values

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('           Reference values\n');
disp('Measure');
fprintf('---------------------------------------------------------------------------------------\n')

ref = getReferenceValues(reference(1:2:end), reference(2:2:end), 'MHR');
fprintf('MHR             '); fprintf('%.2f (%.2f-%.2f)       ',ref); fprintf('\n') 

ref = getReferenceValues(reference(1:2:end), reference(2:2:end), 'SDNN');
fprintf('SDNN            '); fprintf('%.2f (%.2f-%.2f)       ',ref); fprintf('\n')

ref = getReferenceValues(reference(1:2:end), reference(2:2:end), 'RMSSD');
fprintf('RMSSD            '); fprintf('%.2f (%.2f-%.2f)        ',ref); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')

%% Results

errorThreshold = 0.001:0.001:0.25;
% errorThreshold = 0.001:0.001:0.05;
% errorThreshold = 0.01;

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Absolute error: Time indexes, AW\n');
disp('Measure          Method');
fprintf('                 OR           L           NL\n');
fprintf('---------------------------------------------------------------------------------------\n')

figure(1)
error = computeTimeCoverage(results, errorThreshold, 'MHR');
fprintf('MHR              '); fprintf('%.2f (%.2f -- %.2f) &  ',error); fprintf('\n')
% exportgraphics(gca,'aw_mhr.pdf') 

figure(2)
error = computeTimeCoverage(results, errorThreshold, 'SDNN');
fprintf('SDNN             '); fprintf('%.2f (%.2f -- %.2f) &  ',error); fprintf('\n')
% exportgraphics(gca,'aw_sdnn.pdf') 

figure(3)
error = computeTimeCoverage(results, errorThreshold, 'RMSSD');
fprintf('RMSSD            '); fprintf('%.2f (%.2f -- %.2f) &  ',error); fprintf('\n')
% exportgraphics(gca,'aw_rmssd.pdf') 


%% Sympathovagal balance results
 
figure(4)
groupDiscriminationAW(results, 'MHR');
ylabel('MHR [beats/min]')
exportgraphics(gca,'aw_mhr_groups.pdf') 

figure(5)
groupDiscriminationAW(results, 'SDNN');
ylabel('SDNN [ms]')
exportgraphics(gca,'aw_sdnn_groups.pdf') 

figure(6)
groupDiscriminationAW(results, 'RMSSD');
ylabel('RMSSD [ms]');
exportgraphics(gca,'aw_rmssd_groups.pdf') 

%% Statistical differences between methods

metric = 'RMSSD';

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

figure(1)
significance = twoGroupsDegradation(OR', L', 0, metric);
fprintf('pvalues ro-l      '); fprintf('%.2f       ',significance); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')

figure(2)
significance = twoGroupsDegradation(L', NL', 0, metric);
fprintf('pvalues l-nl      '); fprintf('%.2f       ',significance); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')

figure(3)
significance = twoGroupsDegradation(OR', NL', 0, metric);
fprintf('pvalues ro-nl     '); fprintf('%.2f       ',significance); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')


fprintf('---------------------------------------------------------------------------------------\n')

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
            reference{saveIndex} = lpp(referenceSegment,1,false); %#ok<*SAGROW>
            
            % RO
            RO{saveIndex} = lpp(awSegment,1,true);
                      
            % L
            L{saveIndex} = lpp(gapcorrector(awSegment),1,false);
            
            % NL
            NL{saveIndex} = lpp(gapcorrectorNonLinear(awSegment),1,false);
             
                        
            saveIndex = saveIndex+1;
            fprintf('Done\n');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf('Analyzing %s (stress), Segment %d...',files{kk},jj);
            referenceSegment = tRRH7stress(tRRH7stress>segmentBegin(jj) & tRRH7stress<segmentEnd(jj));
            awSegment = tRRAWstress(tRRAWstress>segmentBegin(jj) & tRRAWstress<segmentEnd(jj));
            
            % Reference
            reference{saveIndex} = lpp(referenceSegment,1,false); %#ok<*SAGROW>
            
            % RO
            RO{saveIndex} = lpp(awSegment,1,true);
                      
            % L
            L{saveIndex} = lpp(gapcorrector(awSegment),1,false);
            
            % NL
            NL{saveIndex} = lpp(gapcorrectorNonLinear(awSegment),1,false);
                 

            saveIndex = saveIndex+1;
            fprintf('Done\n');
        end
    end
end

results = [reference; RO; L; NL];



%% Reference values

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('           Reference values\n');
disp('Measure');
fprintf('---------------------------------------------------------------------------------------\n')

ref = getReferenceValues(reference(1:2:end), reference(2:2:end), 'SD1');
fprintf('SD1             '); fprintf('%.2f (%.2f-%.2f)       ',ref); fprintf('\n')

ref = getReferenceValues(reference(1:2:end), reference(2:2:end), 'SD2');
fprintf('SD2            '); fprintf('%.2f (%.2f-%.2f)       ',ref); fprintf('\n')

ref = getReferenceValues(reference(1:2:end), reference(2:2:end), 'Md');
fprintf('Md             '); fprintf('%.2f (%.2f-%.2f)       ',ref); fprintf('\n')

ref = getReferenceValues(reference(1:2:end), reference(2:2:end), 'Sd');
fprintf('Sd            '); fprintf('%.2f (%.2f-%.2f)       ',ref); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')




%% Results

errorThreshold = 0.001:0.001:0.25;
% errorThreshold = 0.01;

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Absolute error: LPP indexes, AW\n');
disp('Measure          Method');
fprintf('                 OR            L            NL\n');
fprintf('---------------------------------------------------------------------------------------\n')

figure(1)
error = computeTimeCoverage(results, errorThreshold, 'SD1');
fprintf('SD1              '); fprintf('%.2f (%.2f -- %.2f)  & ',error); fprintf('\n')
% exportgraphics(gca,'aw_sd1.pdf') 

figure(2)
error = computeTimeCoverage(results, errorThreshold, 'SD2');
fprintf('SD2             '); fprintf('%.2f (%.2f -- %.2f)  & ',error); fprintf('\n')
% exportgraphics(gca,'aw_sd2.pdf') 

figure(3)
error = computeTimeCoverage(results, errorThreshold, 'Md');
fprintf('Md            '); fprintf('%.2f (%.2f -- %.2f)  & ',error); fprintf('\n')
% exportgraphics(gca,'aw_md.pdf') 

figure(4)
error = computeTimeCoverage(results, errorThreshold, 'Sd');
fprintf('Sd            '); fprintf('%.2f (%.2f -- %.2f)  & ',error); fprintf('\n')
% exportgraphics(gca,'aw_sd.pdf') 


%% Sympathovagal balance results
 
figure(1)
groupDiscriminationAW(results, 'SD1');
ylabel('SD1 [ms]')
% exportgraphics(gca,'aw_sd1_groups.pdf') 

figure(2)
groupDiscriminationAW(results, 'SD2');
ylabel('SD2 [ms]')
% exportgraphics(gca,'aw_sd2_groups.pdf') 

figure(3)
groupDiscriminationAW(results, 'SD12');
ylabel('SD12');
% exportgraphics(gca,'aw_sd12_groups.pdf') 

figure(4)
groupDiscriminationAW(results, 'S');
ylabel('S [ms^2]','interpreter','tex');
% exportgraphics(gca,'aw_s_groups.pdf') 

figure(5)
groupDiscriminationAW(results, 'Md');
ylabel('Md [ms]');
% exportgraphics(gca,'aw_md_groups.pdf') 

figure(6)
groupDiscriminationAW(results, 'Sd');
ylabel('Sd [ms]');
% exportgraphics(gca,'aw_sd_groups.pdf') 


%% Statistical differences between methods

metric = 'Sd';

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

figure(1)
significance = twoGroupsDegradation(RO', L', 0, metric);
fprintf('pvalues ro-l      '); fprintf('%.2f       ',significance); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')

figure(2)
significance = twoGroupsDegradation(L', NL', 0, metric);
fprintf('pvalues l-nl      '); fprintf('%.2f       ',significance); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')

figure(3)
significance = twoGroupsDegradation(RO', NL', 0, metric);
fprintf('pvalues ro-nl     '); fprintf('%.2f       ',significance); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')


fprintf('---------------------------------------------------------------------------------------\n')


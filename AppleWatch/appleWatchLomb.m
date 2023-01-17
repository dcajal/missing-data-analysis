%% Lomb's method results with random distributed errors

clear
addpath('lib','AppleWatch/database');
figurePresets

segmentDurationSec = 120;
stepDurationSec = 50;
fs = 4;
nfft = 2^10;
wdw = 60; % seconds
noverlap = 30; % seconds
f = linspace(0,0.4,nfft);
plotflag = false;

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
    segmentEnd = segmentBegin+segmentDurationSec+1;
    
    if ~isempty(segmentBegin)
        for jj = 1:3%length(segmentBegin)  
            fprintf('Analyzing %s (relax), Segment %d...',files{kk},jj);
            referenceSegment = tRRH7relax(tRRH7relax>segmentBegin(jj) & tRRH7relax<segmentEnd(jj));
            awSegment = tRRAWrelax(tRRAWrelax>segmentBegin(jj) & tRRAWrelax<segmentEnd(jj));
                    
            % Reference
            dReference = diff(referenceSegment);
            hrReference = 1000./dReference;
            hrReference = detrend(hrReference);
            psdReference = plombav(referenceSegment(2:end),hrReference,f,wdw,noverlap);
            
                                  
            % Apple Watch L
            tnL = gapcorrector(awSegment);
            dtnL = diff(tnL);
            hrL = 1000./dtnL;
            hrL = detrend(hrL);
            psdL = plombav(tnL(2:end),hrL,f,wdw,noverlap);
            
            % Apple Watch NL
            tnNL = gapcorrectorNonLinear(awSegment);
            dtnNL = diff(tnNL);
            hrNL = 1000./dtnNL;
            hrNL = detrend(hrNL);
            psdNL = plombav(tnNL(2:end),hrNL,f,wdw,noverlap);
            
            % Apple Watch OR
            dAw = diff(awSegment);
            awSegment = awSegment(2:end);
            threshold = computeThreshold(dAw);
            awSegment(dAw>threshold) = [];
            dAw(dAw>threshold) = [];
            hrAw = 1000./dAw;
            hrAw = detrend(hrAw);
            psdOR = plombav(awSegment,hrAw,f,wdw,noverlap);
            
            
            [bb, aa] = butter(15, 0.04*2, 'high');
            h = freqz(bb,aa,nfft);
            psdReference = psdReference.*abs(h);
            psdOR = psdOR.*abs(h);
            psdL = psdL.*abs(h);
            psdNL = psdNL.*abs(h);
            clear h bb aa
            
              
            reference{saveIndex} = freqind(psdReference, f); %#ok<*SAGROW>
            OR{saveIndex} = freqind(psdOR, f); %#ok<*SAGROW>
            L{saveIndex} = freqind(psdL, f); %#ok<*SAGROW>
            NL{saveIndex} = freqind(psdNL, f); %#ok<*SAGROW>
      
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plots

            if plotflag
                figure('DefaultAxesFontSize',14,'units','normalized','outerposition',[0 0 1 1]); 
                subplot(211)
                plot(referenceSegment(2:end),dReference,'Linewidth',2); hold on;
                plot(awSegment,dAw,'Linewidth',1);
                plot(tnL(2:end),dtnL,'Linewidth',1);
                plot(tnNL(2:end),dtnNL,'Linewidth',1);
                legend('Reference (ECG)','OR','L','NL')
                ylabel('NN [s]');
                xlabel('Time [s]');
                subplot(212)
                plot(f,psdReference,'Linewidth',2); hold on;
                plot(f,psdOR,'Linewidth',1);
                plot(f,psdL,'Linewidth',1);
                plot(f,psdNL,'Linewidth',1);
                legend('Reference (ECG)','OR','L','NL');    
                xlabel('Frequency [Hz]')
                ylabel('PSD_{Lomb} [mHz^2]','interpreter','tex')
                axis tight
                pause                
            end
            
            saveIndex = saveIndex+1;
            fprintf('Done\n');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            fprintf('Analyzing %s (stress), Segment %d...',files{kk},jj);
            referenceSegment = tRRH7stress(tRRH7stress>segmentBegin(jj) & tRRH7stress<segmentEnd(jj));
            awSegment = tRRAWstress(tRRAWstress>segmentBegin(jj) & tRRAWstress<segmentEnd(jj));

            % Reference
            dReference = diff(referenceSegment);
            hrReference = 1000./dReference;
            hrReference = detrend(hrReference);
            psdReference = plombav(referenceSegment(2:end),hrReference,f,wdw,noverlap);
            
                                  
            % Apple Watch L
            tnL = gapcorrector(awSegment);
            dtnL = diff(tnL);
            hrL = 1000./dtnL;
            hrL = detrend(hrL);
            psdL = plombav(tnL(2:end),hrL,f,wdw,noverlap);
            
            % Apple Watch NL
            tnNL = gapcorrectorNonLinear(awSegment);
            dtnNL = diff(tnNL);
            hrNL = 1000./dtnNL;
            hrNL = detrend(hrNL);
            psdNL = plombav(tnNL(2:end),hrNL,f,wdw,noverlap);
            
            % Apple Watch OR
            dAw = diff(awSegment);
            awSegment = awSegment(2:end);
            threshold = computeThreshold(dAw);
            awSegment(dAw>threshold) = [];
            dAw(dAw>threshold) = [];
            hrAw = 1000./dAw;
            hrAw = detrend(hrAw);
            psdOR = plombav(awSegment,hrAw,f,wdw,noverlap);
            
              
            reference{saveIndex} = freqind(psdReference, f); %#ok<*SAGROW>
            OR{saveIndex} = freqind(psdOR, f); %#ok<*SAGROW>
            L{saveIndex} = freqind(psdL, f); %#ok<*SAGROW>
            NL{saveIndex} = freqind(psdNL, f); %#ok<*SAGROW>
      
            
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

ref = getReferenceValues(reference(1:2:end), reference(2:2:end), 'LF');
fprintf('LF             '); fprintf('%.2f (%.2f-%.2f)       ',ref); fprintf('\n')

ref = getReferenceValues(reference(1:2:end), reference(2:2:end), 'HF');
fprintf('HF            '); fprintf('%.2f (%.2f-%.2f)       ',ref); fprintf('\n')


fprintf('---------------------------------------------------------------------------------------\n')



%% Results

errorThreshold = 0.001:0.001:0.25;
% errorThreshold = 0.01;

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Absolute error: Freq indexes (Lomb), AW\n');
disp('Measure          Method');
fprintf('                 Remove outliers            Iterative            Iterative NL\n');
fprintf('---------------------------------------------------------------------------------------\n')

figure(1)
error = computeLombCoverage(results, errorThreshold, 'LF');
fprintf('LF              '); fprintf('%.2f (%.2f -- %.2f) &  ',error); fprintf('\n')
% exportgraphics(gca,'aw_lomb_lf.pdf') 

figure(2)
error = computeLombCoverage(results, errorThreshold, 'HF');
fprintf('HF              '); fprintf('%.2f (%.2f -- %.2f) &  ',error); fprintf('\n')
% exportgraphics(gca,'aw_lomb_hf.pdf') 



%% Sympathovagal balance results
 
figure(1)
groupDiscriminationAW(results, 'LF');
ylabel('P_{LF} [mHz^2]','interpreter','tex')
% exportgraphics(gca,'aw_lf_lomb_groups.pdf') 

figure(2)
groupDiscriminationAW(results, 'HF');
ylabel('P_{HF} [mHz^2]','interpreter','tex')
% exportgraphics(gca,'aw_hf_lomb_groups.pdf') 

figure(3)
groupDiscriminationAW(results, 'LFn');
ylabel('P_{LFn}');
% exportgraphics(gca,'aw_lfn_lomb_groups.pdf') 

figure(4)
groupDiscriminationAW(results, 'LFHF');
ylabel('P_{LF}/P_{HF}');
% exportgraphics(gca,'aw_lfhf_lomb_groups.pdf') 


%% Statistical differences between methods

metric = 'HF';

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


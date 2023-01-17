%% Welch's method results with random distributed errors

clear
addpath('lib','AppleWatch/database');
figurePresets

segmentDurationSec = 120;
stepDurationSec = 50;
fs = 4;
nfft = 2^10;
wdw = hamming(60*fs);
noverlap = 30*fs;
f = linspace(0,0.4,nfft);
splineOrder = 4;
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
            [tmReference,idsReference] = incidences(referenceSegment,1);
            tReference = tmReference(1):0.25:tmReference(end);
            ihrReferencesp = mti(tmReference,idsReference,splineOrder);
            ihrReference = spval(ihrReferencesp,tReference);
            ihrReference = ihrReference*1000;
            
            % Apple Watch
            [tmAw,idsAw] = incidences(awSegment,1);
            taw = tmAw(1):0.25:tmAw(end);
            
            ihrMsp = mti(tmAw,idsAw,splineOrder);
            ihrM = spval(ihrMsp,taw);
            ihrM = ihrM*1000;
                    
            ihrLsp = mti(gapcorrector(awSegment),[],splineOrder);
            ihrL = spval(ihrLsp,taw);
            ihrL = ihrL*1000;
            
            ihrNLsp = mti(gapcorrectorNonLinear(awSegment),[],splineOrder);
            ihrNL = spval(ihrNLsp,taw);
            ihrNL = ihrNL*1000;

            
            % Remove baseline
            [bb, aa] = butter(4, 0.04*2/fs, 'high');
            iftReferenceAux= filtfilt(bb, aa, ihrReference);
            ihrMAux = filtfilt(bb, aa, ihrM);
            ihrLAux = filtfilt(bb, aa, ihrL);
            ihrNLAux = filtfilt(bb, aa, ihrNL);

             
            % Compute spectral powers          
            psdReference = pwelch(iftReferenceAux,wdw,noverlap,f,fs);
            reference{saveIndex} = freqind(psdReference,f); %#ok<*SAGROW>

            psdM = pwelch(ihrMAux,wdw,noverlap,f,fs);
            M{saveIndex} = freqind(psdM,f); %#ok<*SAGROW>
          
            psdL = pwelch(ihrLAux,wdw,noverlap,f,fs);
            L{saveIndex} = freqind(psdL,f); %#ok<*SAGROW>
                    
            psdNL = pwelch(ihrNLAux,wdw,noverlap,f,fs);
            NL{saveIndex} = freqind(psdNL,f); %#ok<*SAGROW>
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plots

            if plotflag
                figure('DefaultAxesFontSize',14,'units','normalized','outerposition',[0 0 1 1]); 
                subplot(211)
                plot(tReference,ihrReference,'Linewidth',2); hold on;
                plot(taw,ihrM,'Linewidth',1);
                plot(taw,ihrL,'Linewidth',1);
                plot(taw,ihrNL,'Linewidth',1);
                legend('Reference (ECG)','M','L','NL')
                xlabel('Time [s]')
                ylabel('IHR [mHz]')
                axis tight
                subplot(212)
                plot(f,psdReference,'Linewidth',2); hold on;
                plot(f,psdM,'Linewidth',1);
                plot(f,psdL,'Linewidth',1);
                plot(f,psdNL,'Linewidth',1);
                legend('Reference (ECG)','M','L','NL')    
                xlabel('Frequency [Hz]')
                ylabel('PSD_{Welch} [mHz^2]','interpreter','tex')
                axis tight
                pause                
            end


            saveIndex = saveIndex+1;
            fprintf('Done\n');


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            fprintf('Analyzing %s (stress), Segment %d...',files{kk},jj);
            referenceSegment = tRRH7stress(tRRH7stress>segmentBegin(jj) & tRRH7stress<segmentEnd(jj));
            awSegment = tRRAWstress(tRRAWstress>segmentBegin(jj) & tRRAWstress<segmentEnd(jj));


            % Reference
            [tmReference,idsReference] = incidences(referenceSegment,1);
            tReference = tmReference(1):0.25:tmReference(end);
            ihrReferencesp = mti(tmReference,idsReference,splineOrder);
            ihrReference = spval(ihrReferencesp,tReference);            
            ihrReference = ihrReference*1000;
            
            % Apple Watch
            [tmAw,idsAw,~,tnAw] = incidences(awSegment,1);
            taw = tmAw(1):0.25:tmAw(end);

            
            ihrMsp = mti(tmAw,idsAw,splineOrder);
            ihrM = spval(ihrMsp,taw);
            ihrM = ihrM*1000;
                       
            ihrLsp = mti(gapcorrector(awSegment),[],splineOrder);
            ihrL = spval(ihrLsp,taw);
            ihrL = ihrL*1000;
            
            ihrNLsp = mti(gapcorrectorNonLinear(awSegment),[],splineOrder);
            ihrNL = spval(ihrNLsp,taw);
            ihrNL = ihrNL*1000;
            
            
            % Remove baseline
            [bb, aa] = butter(4, 0.04*2/fs, 'high');
            iftReferenceAux= filtfilt(bb, aa, ihrReference);
            ihrMAux = filtfilt(bb, aa, ihrM);
            ihrLAux = filtfilt(bb, aa, ihrL);
            ihrNLAux = filtfilt(bb, aa, ihrNL);

     
            % Compute spectral powers          
            psdReference = pwelch(iftReferenceAux,wdw,noverlap,f,fs);
            reference{saveIndex} = freqind(psdReference,f); %#ok<*SAGROW>

            psdM = pwelch(ihrMAux,wdw,noverlap,f,fs);
            M{saveIndex} = freqind(psdM,f); %#ok<*SAGROW>
           
            psdL = pwelch(ihrLAux,wdw,noverlap,f,fs);
            L{saveIndex} = freqind(psdL,f); %#ok<*SAGROW>
                    
            psdNL = pwelch(ihrNLAux,wdw,noverlap,f,fs);
            NL{saveIndex} = freqind(psdNL,f); %#ok<*SAGROW>
      
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plots

            if plotflag
                figure('DefaultAxesFontSize',14,'units','normalized','outerposition',[0 0 1 1]); 
                subplot(211)
                plot(tReference,ihrReference,'Linewidth',2); hold on;
                    plot(taw,ihrM,'Linewidth',1);
                    plot(taw,ihrL,'Linewidth',1);
                    plot(taw,ihrNL,'Linewidth',1);
                legend('Reference (ECG)','M','L','NL')
                xlabel('Time [s]')
                ylabel('IHR [mHz]')
                subplot(212)
                plot(f,psdReference,'Linewidth',2); hold on;
                plot(f,psdM,'Linewidth',1);
                plot(f,psdL,'Linewidth',1);
                plot(f,psdNL,'Linewidth',1);
                legend('Reference (ECG)','M','L','NL')    
                xlabel('Frequency [Hz]')
                ylabel('PSD_{Welch} [mHz^2]','interpreter','tex')
                axis tight
                pause                
            end
            
            saveIndex = saveIndex+1;
            fprintf('Done\n');
        end
    end
end

results = [reference; M; L; NL];


%% Reference values

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('           Reference values\n');
disp('Measure');
fprintf('---------------------------------------------------------------------------------------\n')

ref = getReferenceValues(reference(1:2:end), reference(2:2:end), 'LF');
fprintf('P_LF             '); fprintf('%.2f (%.2f-%.2f)       ',ref); fprintf('\n')

ref = getReferenceValues(reference(1:2:end), reference(2:2:end), 'HF');
fprintf('P_HF            '); fprintf('%.2f (%.2f-%.2f)       ',ref); fprintf('\n')


fprintf('---------------------------------------------------------------------------------------\n')

%% Results

errorThreshold = 0.001:0.001:0.25;
% errorThreshold = 0.01;

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Absolute error: Freq indexes (Welch), AW\n');
disp('Measure          Method');
fprintf('                 M            L            NL\n');
fprintf('---------------------------------------------------------------------------------------\n')

figure(1)
error = computeWelchCoverage(results, errorThreshold, 'LF');
% fancyBoxplotAW2(results,'LF')
fprintf('LF              '); fprintf('%.2f (%.2f -- %.2f) &  ',error); fprintf('\n')
% exportgraphics(gca,'aw_lf.pdf') 

figure(2)
error = computeWelchCoverage(results, errorThreshold, 'HF');
% fancyBoxplotAW2(results,'HF')
fprintf('HF              '); fprintf('%.2f (%.2f -- %.2f) &  ',error); fprintf('\n')
% exportgraphics(gca,'aw_hf.pdf')


%% Sympathovagal balance results
 
figure(1)
groupDiscriminationAW(results, 'LF');
ylabel('P_{LF} [mHz^2]','interpreter','tex')
% exportgraphics(gca,'aw_lf_groups.pdf') 

figure(2)
groupDiscriminationAW(results, 'HF');
ylabel('P_{HF} [mHz^2]','interpreter','tex')
% exportgraphics(gca,'aw_hf_groups.pdf') 

figure(3)
groupDiscriminationAW(results, 'LFn');
ylabel('LFn');
% exportgraphics(gca,'aw_lfn_groups.pdf') 

figure(4)
groupDiscriminationAW(results, 'LFHF');
ylabel('LF/HF');
% exportgraphics(gca,'aw_lfhf_groups.pdf') 


%% Statistical differences between methods

metric = 'HF';

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

figure(1)
significance = twoGroupsDegradation(M', L', 0, metric);
fprintf('pvalues m-l       '); fprintf('%.2f       ',significance); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')

figure(2)
significance = twoGroupsDegradation(L', NL', 0, metric);
fprintf('pvalues l-nl      '); fprintf('%.2f       ',significance); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')

figure(3)
significance = twoGroupsDegradation(M', NL', 0, metric);
fprintf('pvalues m-nl      '); fprintf('%.2f       ',significance); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')


fprintf('---------------------------------------------------------------------------------------\n')

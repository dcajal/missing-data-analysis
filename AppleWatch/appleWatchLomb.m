%% Lomb's method results with random distributed errors

clear
addpath('lib','AppleWatch/database');
figurePresets

segmentDurationSec = 2*60;
stepDurationSec = 1*60;
fs = 4;
nfft = 2^10;
wdw = 60; % seconds
noverlap = 30; % seconds
f = linspace(0,0.4,nfft);

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
            dReferenceSegment = diff(referenceSegment);
            dReferenceSegment = 1000.*dReferenceSegment;
            dReferenceSegment = detrend(dReferenceSegment);
            psdReference = plombav(referenceSegment(2:end),dReferenceSegment,f,wdw,noverlap);
            
            % Apple Watch none
            dAw = diff(awSegment);
            threshold = computeThreshold(dAw);
            awSegment(dAw>threshold) = [];
            dAw(dAw>threshold) = [];
            dAw = 1000.*dAw;
            dAw = detrend(dAw);
            psdAwNone = plombav(awSegment(2:end),dAw,f,wdw,noverlap);
            
            % Apple Watch incidences
            [~,~,~,tnAwIncidences] = incidences(awSegment);
            dtnAwIncidences = diff(tnAwIncidences);
            dtnAwIncidences = 1000.*dtnAwIncidences;
            dtnAwIncidences = detrend(dtnAwIncidences);
            psdAwIncidences = plombav(tnAwIncidences(2:end),dtnAwIncidences,f,wdw,noverlap);
            
            % Apple Watch iterative
            tnAwIterative = gapcorrector(awSegment);
            dtnAwIterative = diff(tnAwIterative);
            dtnAwIterative = 1000.*dtnAwIterative;
            dtnAwIterative = detrend(dtnAwIterative);
            psdAwIterative = plombav(tnAwIterative(2:end),dtnAwIterative,f,wdw,noverlap);
            
            % Apple Watch iterative non linear
            tnAwIterativeNL = gapcorrectorNonLinear(awSegment);
            dtnAwIterativeNL = diff(tnAwIterativeNL);
            dtnAwIterativeNL = 1000.*dtnAwIterativeNL;
            dtnAwIterativeNL = detrend(dtnAwIterativeNL);
            psdAwIterativeNL = plombav(tnAwIterativeNL(2:end),dtnAwIterativeNL,f,wdw,noverlap);
            
            [bb, aa] = butter(15, 0.04*2, 'high');
            h = freqz(bb,aa,nfft);
            psdReference = psdReference.*abs(h);
            psdAwNone = psdAwNone.*abs(h);
            psdAwIncidences = psdAwIncidences.*abs(h);
            psdAwIterative = psdAwIterative.*abs(h);
            psdAwIterativeNL = psdAwIterativeNL.*abs(h);
            clear h bb aa
            
              
            reference{saveIndex} = freqind(psdReference, f); %#ok<*SAGROW>
            awNone{saveIndex} = freqind(psdAwNone, f); %#ok<*SAGROW>
            awIncidences{saveIndex} = freqind(psdAwIncidences, f); %#ok<*SAGROW>
            awIterative{saveIndex} = freqind(psdAwIterative, f); %#ok<*SAGROW>
            awIterativeNL{saveIndex} = freqind(psdAwIterativeNL, f); %#ok<*SAGROW>
      
            
            saveIndex = saveIndex+1;
            fprintf('Done\n');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            fprintf('Analyzing %s (stress), Segment %d...',files{kk},jj);
            referenceSegment = tRRH7stress(tRRH7stress>segmentBegin(jj) & tRRH7stress<segmentEnd(jj));
            awSegment = tRRAWstress(tRRAWstress>segmentBegin(jj) & tRRAWstress<segmentEnd(jj));

            % Reference
            dReferenceSegment = diff(referenceSegment);
            dReferenceSegment = 1000.*dReferenceSegment;
            dReferenceSegment = detrend(dReferenceSegment);
            psdReference = plombav(referenceSegment(2:end),dReferenceSegment,f,wdw,noverlap);
            
            % Apple Watch none
            dAw = diff(awSegment);
            threshold = computeThreshold(dAw);
            awSegment(dAw>threshold) = [];
            dAw(dAw>threshold) = [];
            dAw = 1000.*dAw;
            dAw = detrend(dAw);
            psdAwNone = plombav(awSegment(2:end),dAw,f,wdw,noverlap);
            
            % Apple Watch incidences
            [~,~,~,tnAwIncidences] = incidences(awSegment);
            dtnAwIncidences = diff(tnAwIncidences);
            dtnAwIncidences = 1000.*dtnAwIncidences;
            dtnAwIncidences = detrend(dtnAwIncidences);
            psdAwIncidences = plombav(tnAwIncidences(2:end),dtnAwIncidences,f,wdw,noverlap);
            
            % Apple Watch iterative
            tnAwIterative = gapcorrector(awSegment);
            dtnAwIterative = diff(tnAwIterative);
            dtnAwIterative = 1000.*dtnAwIterative;
            dtnAwIterative = detrend(dtnAwIterative);
            psdAwIterative = plombav(tnAwIterative(2:end),dtnAwIterative,f,wdw,noverlap);
            
            % Apple Watch iterative non linear
            tnAwIterativeNL = gapcorrectorNonLinear(awSegment);
            dtnAwIterativeNL = diff(tnAwIterativeNL);
            dtnAwIterativeNL = 1000.*dtnAwIterativeNL;
            dtnAwIterativeNL = detrend(dtnAwIterativeNL);
            psdAwIterativeNL = plombav(tnAwIterativeNL(2:end),dtnAwIterativeNL,f,wdw,noverlap);
            
            [bb, aa] = butter(15, 0.04*2, 'high');
            h = freqz(bb,aa,nfft);
            psdReference = psdReference.*abs(h);
            psdAwNone = psdAwNone.*abs(h);
            psdAwIncidences = psdAwIncidences.*abs(h);
            psdAwIterative = psdAwIterative.*abs(h);
            psdAwIterativeNL = psdAwIterativeNL.*abs(h);
            clear h bb aa
            
              
            reference{saveIndex} = freqind(psdReference, f); %#ok<*SAGROW>
            awNone{saveIndex} = freqind(psdAwNone, f); %#ok<*SAGROW>
            awIncidences{saveIndex} = freqind(psdAwIncidences, f); %#ok<*SAGROW>
            awIterative{saveIndex} = freqind(psdAwIterative, f); %#ok<*SAGROW>
            awIterativeNL{saveIndex} = freqind(psdAwIterativeNL, f); %#ok<*SAGROW>
      
            
            saveIndex = saveIndex+1;
            fprintf('Done\n');
        end
    end
end

results = [reference; awNone; awIncidences; awIterative; awIterativeNL];


%             figure;
%             plot(tReference,mtReference,'-'); hold on;
%             plot(taw,mtAwNone,'--');
%             plot(taw,mtAwIncidences,'-*');
%             plot(taw,mtAwIterative,'-^');
%             legend('Reference','IPFM','Incidences','Iterative')
%             axis tight; set(gcf, 'Position', get(0, 'Screensize'));
%             pause;
            
%             figure;
%             plot(f,psdReference,'-'); hold on;
%             plot(f,psdAwIpfm,'--');
%             plot(f,psdAwIncidences,'-*');
%             plot(f,psdAwIterative,'-^');
%             legend('Reference','IPFM','Incidences','Iterative')
%             axis tight; set(gcf, 'Position', get(0, 'Screensize'));
%             pause;

%% Results

errorThreshold = 0.001:0.001:0.3;
% errorThreshold = 0.01;

computeLombCoverage(results, errorThreshold, 'HF');
computeLombCoverage(results, errorThreshold, 'LF');
computeLombCoverage(results, errorThreshold, 'LFn');
computeLombCoverage(results, errorThreshold, 'LFHF');

%% Welch's method results with random distributed errors

clear
addpath('lib','AppleWatch/database');
figurePresets

segmentDurationSec = 2*60;
stepDurationSec = 1*60;
fs = 4;
nfft = 2^10;
wdw = hamming(60*fs);
noverlap = 30*fs;
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
            
            tReference = referenceSegment(1):0.25:referenceSegment(end);
            taw = awSegment(1):0.25:awSegment(end);
            
            
            % Reference (incidences)
            [tmReference,idsReference] = incidences(referenceSegment);
            mtReference = ipfm(tmReference,idsReference,tReference);
            
            % Apple Watch
            [tmAw,idsAw,~,tnAw] = incidences(awSegment);
            mtAwIpfm = ipfm(tmAw,idsAw,taw);
            mtAwIncidences = ipfm(tnAw,[],taw);
            mtAwIterative = ipfm(gapcorrector(awSegment),[],taw);
            mtAwIterativeNL = ipfm(gapcorrectorNonLinear(awSegment),[],taw);

            
            % Remove baseline
            [bb, aa] = butter(4, 0.04*2/fs, 'high');
            auxMtReference= filtfilt(bb, aa, 1000./mtReference);
            auxMtAwIpfm = filtfilt(bb, aa, 1000./mtAwIpfm);
            auxMtAwIncidences = filtfilt(bb, aa, 1000./mtAwIncidences);
            auxMtAwIterative = filtfilt(bb, aa, 1000./mtAwIterative);
            auxMtAwIterativeNL = filtfilt(bb, aa, 1000./mtAwIterativeNL);

             
            % Compute spectral powers          
            psdReference = pwelch(auxMtReference,wdw,noverlap,f,fs);
            reference{saveIndex} = freqind(psdReference,f); %#ok<*SAGROW>

            psdAwIpfm = pwelch(auxMtAwIpfm,wdw,noverlap,f,fs);
            awIpfm{saveIndex} = freqind(psdAwIpfm,f); %#ok<*SAGROW>
           
            psdAwIncidences = pwelch(auxMtAwIncidences,wdw,noverlap,f,fs);
            awIncidences{saveIndex} = freqind(psdAwIncidences,f); %#ok<*SAGROW>

            psdAwIterative = pwelch(auxMtAwIterative,wdw,noverlap,f,fs);
            awIterative{saveIndex} = freqind(psdAwIterative,f); %#ok<*SAGROW>
                    
            psdAwIterativeNL = pwelch(auxMtAwIterativeNL,wdw,noverlap,f,fs);
            awIterativeNL{saveIndex} = freqind(psdAwIterativeNL,f); %#ok<*SAGROW>


            saveIndex = saveIndex+1;
            fprintf('Done\n');


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            fprintf('Analyzing %s (stress), Segment %d...',files{kk},jj);
            referenceSegment = tRRH7stress(tRRH7stress>segmentBegin(jj) & tRRH7stress<segmentEnd(jj));
            awSegment = tRRAWstress(tRRAWstress>segmentBegin(jj) & tRRAWstress<segmentEnd(jj));

            tReference = referenceSegment(1):0.25:referenceSegment(end);
            taw = awSegment(1):0.25:awSegment(end);


            % Reference (incidences)
            [tmReference,idsReference] = incidences(referenceSegment);
            mtReference = ipfm(tmReference,idsReference,tReference);
            
            % Apple Watch
            [tmAw,idsAw,~,tnAw] = incidences(awSegment);
            mtAwIpfm = ipfm(tmAw,idsAw,taw);
            mtAwIncidences = ipfm(tnAw,[],taw);
            mtAwIterative = ipfm(gapcorrector(awSegment),[],taw);
            mtAwIterativeNL = ipfm(gapcorrectorNonLinear(awSegment),[],taw);
            
            
            % Remove baseline
            [bb, aa] = butter(4, 0.04*2/fs, 'high');
            auxMtReference= filtfilt(bb, aa, 1000./mtReference);
            auxMtAwIpfm = filtfilt(bb, aa, 1000./mtAwIpfm);
            auxMtAwIncidences = filtfilt(bb, aa, 1000./mtAwIncidences);
            auxMtAwIterative = filtfilt(bb, aa, 1000./mtAwIterative);
            auxMtAwIterativeNL = filtfilt(bb, aa, 1000./mtAwIterativeNL);

     
            % Compute spectral powers          
            psdReference = pwelch(auxMtReference,wdw,noverlap,f,fs);
            reference{saveIndex} = freqind(psdReference,f); %#ok<*SAGROW>

            psdAwIpfm = pwelch(auxMtAwIpfm,wdw,noverlap,f,fs);
            awIpfm{saveIndex} = freqind(psdAwIpfm,f); %#ok<*SAGROW>
           
            psdAwIncidences = pwelch(auxMtAwIncidences,wdw,noverlap,f,fs);
            awIncidences{saveIndex} = freqind(psdAwIncidences,f); %#ok<*SAGROW>

            psdAwIterative = pwelch(auxMtAwIterative,wdw,noverlap,f,fs);
            awIterative{saveIndex} = freqind(psdAwIterative,f); %#ok<*SAGROW>
                    
            psdAwIterativeNL = pwelch(auxMtAwIterativeNL,wdw,noverlap,f,fs);
            awIterativeNL{saveIndex} = freqind(psdAwIterativeNL,f); %#ok<*SAGROW>
      
           
            saveIndex = saveIndex+1;
            fprintf('Done\n');
        end
    end
end

results = [reference; awIpfm; awIncidences; awIterative; awIterativeNL];


%             figure;
%             plot(tReference,mtReference,'-'); hold on;
%             plot(taw,mtAwIpfm,'--');
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

computeWelchCoverage(results, errorThreshold, 'HF');
computeWelchCoverage(results, errorThreshold, 'LF');
computeWelchCoverage(results, errorThreshold, 'LFn');
computeWelchCoverage(results, errorThreshold, 'LFHF');


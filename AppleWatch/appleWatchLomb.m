%% Lomb's method results with random distributed errors

clear
addpath('lib','database');
set(0,'defaultAxesFontName', 'Helvetica');
set(0,'defaultTextFontName', 'Helvetica');
set(0,'defaulttextInterpreter','latex')
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize', 20);
set(0,'defaultTextFontSize', 20);
set(0,'defaultFigureRenderer', 'painters')

segmentDurationSec = 2*60;
stepDurationSec = 1*60;
fs = 4;
nfft = 2^10;
wdw = 60; % seconds
noverlap = 30; % seconds
f = linspace(0,0.4,nfft);

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
        for jj = 1:length(segmentBegin)  
            fprintf('Analyzing %s (relax), Segment %d...',files{kk},jj);
            referenceSegment = tRRH7relax(tRRH7relax>segmentBegin(jj) & tRRH7relax<segmentEnd(jj));
            awSegment = tRRAWrelax(tRRAWrelax>segmentBegin(jj) & tRRAWrelax<segmentEnd(jj));
            
            tReference = referenceSegment(1):0.25:referenceSegment(end);
            taw = awSegment(1):0.25:awSegment(end);
            
            
            % Reference
            dReferenceSegment = diff(referenceSegment);
            dReferenceSegment = 1000.*dReferenceSegment;
            dReferenceSegment = detrend(dReferenceSegment);
            psdReference = plombav(referenceSegment(2:end),dReferenceSegment,f,wdw,noverlap);
            
            
            % Apple Watch incidences
            [~,~,~,tnAwIncidences] = incidences(awSegment);
            dtnAwIncidences = diff(tnAwIncidences);
            dtnAwIncidences = 1000.*dtnAwIncidences;
            dtnAwIncidences = detrend(dtnAwIncidences);
            psdAwIncidences = plombav(tnAwIncidences(2:end),dtnAwIncidences,f,wdw,noverlap);
            
            % Apple Watch iterative
            tnAwIterative = gapcorrectorNonLinear(awSegment);
            dtnAwIterative = diff(tnAwIterative);
            dtnAwIterative = 1000.*dtnAwIterative;
            dtnAwIterative = detrend(dtnAwIterative);
            psdAwIterative = plombav(tnAwIterative(2:end),dtnAwIterative,f,wdw,noverlap);
            
            % Apple Watch none
            dAw = diff(awSegment);
            threshold = computeThreshold(dAw);
            awSegment(dAw>threshold) = [];
            dAw(dAw>threshold) = [];
            dAw = 1000.*dAw;
            dAw = detrend(dAw);
            psdAwNone = plombav(awSegment(2:end),dAw,f,wdw,noverlap);
            
            
            [bb, aa] = butter(15, 0.04*2, 'high');
            h = freqz(bb,aa,nfft);
            psdReference = psdReference.*abs(h);
            psdAwIncidences = psdAwIncidences.*abs(h);
            psdAwIterative = psdAwIterative.*abs(h);
            psdAwNone = psdAwNone.*abs(h);
            clear h bb aa
%             
%             figure;
%             plot(f,psdReference,'-'); hold on;
%             plot(f,psdNone,'--');
%             plot(f,psdIncidences,'-*');
%             plot(f,psdIterative,'-^');
%             legend('Reference','None','Incidences','Iterative')
%             axis tight; set(gcf, 'Position', get(0, 'Screensize'));
%             pause;

            
     
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % HF
            HFBand = [find(f>=0.15,1) find(f>=0.4,1)];
            LFBand = [find(f>=0.04,1) find(f>=0.15,1)];
            
          
            reference{saveIndex}.hf = trapz(f(HFBand(1):HFBand(2)),psdReference(HFBand(1):HFBand(2))); %#ok<*SAGROW>
            reference{saveIndex}.lf = trapz(f(LFBand(1):LFBand(2)),psdReference(LFBand(1):LFBand(2)));
            reference{saveIndex}.lfn = reference{saveIndex}.lf/(reference{saveIndex}.lf+reference{saveIndex}.hf);
            reference{saveIndex}.lfhf = reference{saveIndex}.lf/reference{saveIndex}.hf;
         
            awNone{saveIndex}.hf = trapz(f(HFBand(1):HFBand(2)),psdAwNone(HFBand(1):HFBand(2)));
            awNone{saveIndex}.lf = trapz(f(LFBand(1):LFBand(2)),psdAwNone(LFBand(1):LFBand(2)));
            awNone{saveIndex}.lfn = awNone{saveIndex}.lf/(awNone{saveIndex}.lf+awNone{saveIndex}.hf);
            awNone{saveIndex}.lfhf = awNone{saveIndex}.lf/awNone{saveIndex}.hf;

            awIncidences{saveIndex}.hf = trapz(f(HFBand(1):HFBand(2)),psdAwIncidences(HFBand(1):HFBand(2)));
            awIncidences{saveIndex}.lf = trapz(f(LFBand(1):LFBand(2)),psdAwIncidences(LFBand(1):LFBand(2)));
            awIncidences{saveIndex}.lfn = awIncidences{saveIndex}.lf/(awIncidences{saveIndex}.lf+awIncidences{saveIndex}.hf);
            awIncidences{saveIndex}.lfhf = awIncidences{saveIndex}.lf/awIncidences{saveIndex}.hf;

            awIterative{saveIndex}.hf = trapz(f(HFBand(1):HFBand(2)),psdAwIterative(HFBand(1):HFBand(2)));
            awIterative{saveIndex}.lf = trapz(f(LFBand(1):LFBand(2)),psdAwIterative(LFBand(1):LFBand(2)));
            awIterative{saveIndex}.lfn = awIterative{saveIndex}.lf/(awIterative{saveIndex}.lf+awIterative{saveIndex}.hf);
            awIterative{saveIndex}.lfhf = awIterative{saveIndex}.lf/awIterative{saveIndex}.hf;
         
            
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
        for jj = 1:length(segmentBegin)  
            fprintf('Analyzing %s (stress), Segment %d...',files{kk},jj);
            referenceSegment = tRRH7stress(tRRH7stress>segmentBegin(jj) & tRRH7stress<segmentEnd(jj));
            awSegment = tRRAWstress(tRRAWstress>segmentBegin(jj) & tRRAWstress<segmentEnd(jj));

            tReference = referenceSegment(1):0.25:referenceSegment(end);
            taw = awSegment(1):0.25:awSegment(end);


            % Reference (incidences)
            [tmReference,idsReference] = incidences(referenceSegment);
            mtReference = ipfm(tmReference,idsReference,tReference);

            % Apple Watch
            [tmarmband,idsarmband,~,tnAwIncidences] = incidences(awSegment);
            mtAwNone = ipfm(tmarmband,idsarmband,taw);
            mtAwIncidences = ipfm(tnAwIncidences,[],taw);
            mtAwIterative = ipfm(gapcorrector(awSegment),[],taw);

            
%             figure;
%             plot(tReference,mtReference,'-'); hold on;
%             plot(taw,mtAwNone,'--');
%             plot(taw,mtAwIncidences,'-*');
%             plot(taw,mtAwIterative,'-^');
%             legend('Reference','IPFM','Incidences','Iterative')
%             axis tight; set(gcf, 'Position', get(0, 'Screensize'));
%             pause;
            
            
            % Remove baseline
            [bb, aa] = butter(4, 0.04*2/fs, 'high');
            auxMtReference= filtfilt(bb, aa, 1000./mtReference);
            auxmtAwNone = filtfilt(bb, aa, 1000./mtAwNone);
            auxMtAwIncidences = filtfilt(bb, aa, 1000./mtAwIncidences);
            auxMtAwIterative = filtfilt(bb, aa, 1000./mtAwIterative);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % HF
            HFBand = [find(f>=0.15,1) find(f>=0.4,1)];
            LFBand = [find(f>=0.04,1) find(f>=0.15,1)];

            if ~isempty(auxMtReference)
                psdReference = pwelch(auxMtReference,wdw,noverlap,f,fs);
                reference{saveIndex}.hf = trapz(f(HFBand(1):HFBand(2)),psdReference(HFBand(1):HFBand(2))); %#ok<*SAGROW>
                reference{saveIndex}.lf = trapz(f(LFBand(1):LFBand(2)),psdReference(LFBand(1):LFBand(2)));
                reference{saveIndex}.lfn = reference{saveIndex}.lf/(reference{saveIndex}.lf+reference{saveIndex}.hf);
                reference{saveIndex}.lfhf = reference{saveIndex}.lf/reference{saveIndex}.hf;
            else
                reference{saveIndex}.hf = nan;
                reference{saveIndex}.lf = nan;
                reference{saveIndex}.lfn = nan;
                reference{saveIndex}.lfhf = nan;
            end

            if ~isempty(auxmtAwNone)
                psdAwIpfm = pwelch(auxmtAwNone,wdw,noverlap,f,fs);
                awNone{saveIndex}.hf = trapz(f(HFBand(1):HFBand(2)),psdAwIpfm(HFBand(1):HFBand(2)));
                awNone{saveIndex}.lf = trapz(f(LFBand(1):LFBand(2)),psdAwIpfm(LFBand(1):LFBand(2)));
                awNone{saveIndex}.lfn = awNone{saveIndex}.lf/(awNone{saveIndex}.lf+awNone{saveIndex}.hf);
                awNone{saveIndex}.lfhf = awNone{saveIndex}.lf/awNone{saveIndex}.hf;
            else
                awNone{saveIndex}.hf = nan;
                awNone{saveIndex}.lf = nan;
                awNone{saveIndex}.lfn = nan;
                awNone{saveIndex}.lfhf = nan;
            end

            if ~isempty(auxMtAwIncidences)
                psdAwIncidences = pwelch(auxMtAwIncidences,wdw,noverlap,f,fs);
                awIncidences{saveIndex}.hf = trapz(f(HFBand(1):HFBand(2)),psdAwIncidences(HFBand(1):HFBand(2)));
                awIncidences{saveIndex}.lf = trapz(f(LFBand(1):LFBand(2)),psdAwIncidences(LFBand(1):LFBand(2)));
                awIncidences{saveIndex}.lfn = awIncidences{saveIndex}.lf/(awIncidences{saveIndex}.lf+awIncidences{saveIndex}.hf);
                awIncidences{saveIndex}.lfhf = awIncidences{saveIndex}.lf/awIncidences{saveIndex}.hf;
            else
                awIncidences{saveIndex}.hf = nan;
                awIncidences{saveIndex}.lf = nan;
                awIncidences{saveIndex}.lfn = nan;
                awIncidences{saveIndex}.lfhf = nan;
            end

            if ~isempty(auxMtAwIterative)         
                psdAwIterative = pwelch(auxMtAwIterative,wdw,noverlap,f,fs);
                awIterative{saveIndex}.hf = trapz(f(HFBand(1):HFBand(2)),psdAwIterative(HFBand(1):HFBand(2)));
                awIterative{saveIndex}.lf = trapz(f(LFBand(1):LFBand(2)),psdAwIterative(LFBand(1):LFBand(2)));
                awIterative{saveIndex}.lfn = awIterative{saveIndex}.lf/(awIterative{saveIndex}.lf+awIterative{saveIndex}.hf);
                awIterative{saveIndex}.lfhf = awIterative{saveIndex}.lf/awIterative{saveIndex}.hf;
            else
                awIterative{saveIndex}.hf = nan;
                awIterative{saveIndex}.lf = nan;
                awIterative{saveIndex}.lfn = nan;
                awIterative{saveIndex}.lfhf = nan;
            end
            
%             figure;
%             plot(f,psdReference,'-'); hold on;
%             plot(f,psdAwIpfm,'--');
%             plot(f,psdAwIncidences,'-*');
%             plot(f,psdAwIterative,'-^');
%             legend('Reference','IPFM','Incidences','Iterative')
%             axis tight; set(gcf, 'Position', get(0, 'Screensize'));
%             pause;

            saveIndex = saveIndex+1;
            fprintf('Done\n');
        end
    end
end

%%
errorThreshold = 0.001:0.001:0.1;

%% HF

ipfmCorrects = zeros(length(reference),numel(errorThreshold));
incidencesCorrects = zeros(length(reference),numel(errorThreshold));
iterativeCorrects = zeros(length(reference),numel(errorThreshold));
for kk=1:length(reference)
    ipfmCorrects(kk,:) = abs((awNone{kk}.hf - reference{kk}.hf)/reference{kk}.hf)<errorThreshold;
    incidencesCorrects(kk,:) = abs((awIncidences{kk}.hf - reference{kk}.hf)/reference{kk}.hf)<errorThreshold;
    iterativeCorrects(kk,:) = abs((awIterative{kk}.hf - reference{kk}.hf)/reference{kk}.hf)<errorThreshold;
end

% fprintf('\n');
% fprintf('----------------------------------------\n');
% fprintf('----------------- HF -------------------\n');
% fprintf('----------------------------------------\n');
% fprintf('IPFM coverage: %i out of %i (%.2f%%)\n',sum(ipfmCorrects),numel(ipfmCorrects),sum(ipfmCorrects)/numel(ipfmCorrects)*100);
% fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
% fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
% fprintf('----------------------------------------\n');
% fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
% fprintf('Iterative coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(iterativeCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
% fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
% fprintf('Incidences coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(incidencesCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
% fprintf('----------------------------------------\n');

% figure; hold on
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 3 3],[1 1 1]);
% p(1) = bar(1:length(reference),ipfmCorrects*3,'b');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 2 2],[1 1 1]);
% p(2) = bar(1:length(reference),incidencesCorrects*2,'g');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 1 1],[1 1 1]);
% p(3) = bar(1:length(reference),iterativeCorrects,'k');
% legend(p,'IPFM','Incidences','Iterative','Location','bestoutside')
% title('HF'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

figure; hold on
plot(100*errorThreshold,sum(ipfmCorrects,1),'-');
plot(100*errorThreshold,sum(incidencesCorrects,1),'-*');
plot(100*errorThreshold,sum(iterativeCorrects,1),'-^');
legend('IPFM','Incidences','Iterative','Location','best')
axis tight; set(gcf, 'Position', get(0, 'Screensize'));
title('HF'); xlabel('Permited error (\%)'); ylabel('Number of correct cases')


%% LF

ipfmCorrects = zeros(length(reference),numel(errorThreshold));
incidencesCorrects = zeros(length(reference),numel(errorThreshold));
iterativeCorrects = zeros(length(reference),numel(errorThreshold));
for kk=1:length(reference)
    ipfmCorrects(kk,:) = abs((awNone{kk}.lf - reference{kk}.lf)/reference{kk}.lf)<errorThreshold;
    incidencesCorrects(kk,:) = abs((awIncidences{kk}.lf - reference{kk}.lf)/reference{kk}.lf)<errorThreshold;
    iterativeCorrects(kk,:) = abs((awIterative{kk}.lf - reference{kk}.lf)/reference{kk}.lf)<errorThreshold;
end

% fprintf('\n');
% fprintf('----------------------------------------\n');
% fprintf('----------------- LF -------------------\n');
% fprintf('----------------------------------------\n');
% fprintf('IPFM coverage: %i out of %i (%.2f%%)\n',sum(ipfmCorrects),numel(ipfmCorrects),sum(ipfmCorrects)/numel(ipfmCorrects)*100);
% fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
% fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
% fprintf('----------------------------------------\n');
% fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
% fprintf('Iterative coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(iterativeCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
% fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
% fprintf('Incidences coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(incidencesCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
% fprintf('----------------------------------------\n');

% figure; hold on
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 3 3],[1 1 1]);
% p(1) = bar(1:length(reference),ipfmCorrects*3,'b');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 2 2],[1 1 1]);
% p(2) = bar(1:length(reference),incidencesCorrects*2,'g');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 1 1],[1 1 1]);
% p(3) = bar(1:length(reference),iterativeCorrects,'k');
% legend(p,'IPFM','Incidences','Iterative','Location','bestoutside')
% title('LF'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

figure; hold on
plot(100*errorThreshold,sum(ipfmCorrects,1),'-');
plot(100*errorThreshold,sum(incidencesCorrects,1),'-*');
plot(100*errorThreshold,sum(iterativeCorrects,1),'-^');
legend('IPFM','Incidences','Iterative','Location','best')
axis tight; set(gcf, 'Position', get(0, 'Screensize'));
title('LF'); xlabel('Permited error (\%)'); ylabel('Number of correct cases')


%% LFn

ipfmCorrects = zeros(length(reference),numel(errorThreshold));
incidencesCorrects = zeros(length(reference),numel(errorThreshold));
iterativeCorrects = zeros(length(reference),numel(errorThreshold));
for kk=1:length(reference)
    ipfmCorrects(kk,:) = abs((awNone{kk}.lfn - reference{kk}.lfn)/reference{kk}.lfn)<errorThreshold;
    incidencesCorrects(kk,:) = abs((awIncidences{kk}.lfn - reference{kk}.lfn)/reference{kk}.lfn)<errorThreshold;
    iterativeCorrects(kk,:) = abs((awIterative{kk}.lfn - reference{kk}.lfn)/reference{kk}.lfn)<errorThreshold;
end

% fprintf('\n');
% fprintf('----------------------------------------\n');
% fprintf('---------------- LFn -------------------\n');
% fprintf('----------------------------------------\n');
% fprintf('IPFM coverage: %i out of %i (%.2f%%)\n',sum(ipfmCorrects),numel(ipfmCorrects),sum(ipfmCorrects)/numel(ipfmCorrects)*100);
% fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
% fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
% fprintf('----------------------------------------\n');
% fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
% fprintf('Iterative coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(iterativeCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
% fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
% fprintf('Incidences coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(incidencesCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
% fprintf('----------------------------------------\n');

% figure; hold on
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 3 3],[1 1 1]);
% p(1) = bar(1:length(reference),ipfmCorrects*3,'b');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 2 2],[1 1 1]);
% p(2) = bar(1:length(reference),incidencesCorrects*2,'g');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 1 1],[1 1 1]);
% p(3) = bar(1:length(reference),iterativeCorrects,'k');
% legend(p,'IPFM','Incidences','Iterative','Location','bestoutside')
% title('LFn'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

figure; hold on
plot(100*errorThreshold,sum(ipfmCorrects,1),'-');
plot(100*errorThreshold,sum(incidencesCorrects,1),'-*');
plot(100*errorThreshold,sum(iterativeCorrects,1),'-^');
legend('IPFM','Incidences','Iterative','Location','best')
axis tight; set(gcf, 'Position', get(0, 'Screensize'));
title('LFn'); xlabel('Permited error (\%)'); ylabel('Number of correct cases')



%% LF/HF

ipfmCorrects = zeros(length(reference),numel(errorThreshold));
incidencesCorrects = zeros(length(reference),numel(errorThreshold));
iterativeCorrects = zeros(length(reference),numel(errorThreshold));
for kk=1:length(reference)
    ipfmCorrects(kk,:) = abs((awNone{kk}.lfhf - reference{kk}.lfhf)/reference{kk}.lfhf)<errorThreshold;
    incidencesCorrects(kk,:) = abs((awIncidences{kk}.lfhf - reference{kk}.lfhf)/reference{kk}.lfhf)<errorThreshold;
    iterativeCorrects(kk,:) = abs((awIterative{kk}.lfhf - reference{kk}.lfhf)/reference{kk}.lfhf)<errorThreshold;
end

% fprintf('\n');
% fprintf('----------------------------------------\n');
% fprintf('--------------- LF/HF ------------------\n');
% fprintf('----------------------------------------\n');
% fprintf('IPFM coverage: %i out of %i (%.2f%%)\n',sum(ipfmCorrects),numel(ipfmCorrects),sum(ipfmCorrects)/numel(ipfmCorrects)*100);
% fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
% fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
% fprintf('----------------------------------------\n');
% fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
% fprintf('Iterative coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(iterativeCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
% fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
% fprintf('Incidences coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(incidencesCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
% fprintf('----------------------------------------\n');

% figure; hold on
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 3 3],[1 1 1]);
% p(1) = bar(1:length(reference),ipfmCorrects*3,'b');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 2 2],[1 1 1]);
% p(2) = bar(1:length(reference),incidencesCorrects*2,'g');
% patch([0 length(reference)+1 length(reference)+1 0],[0 0 1 1],[1 1 1]);
% p(3) = bar(1:length(reference),iterativeCorrects,'k');
% legend(p,'IPFM','Incidences','Iterative','Location','bestoutside')
% title('LF/HF'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

figure; hold on
plot(100*errorThreshold,sum(ipfmCorrects,1),'-');
plot(100*errorThreshold,sum(incidencesCorrects,1),'-*');
plot(100*errorThreshold,sum(iterativeCorrects,1),'-^');
legend('IPFM','Incidences','Iterative','Location','best')
axis tight; set(gcf, 'Position', get(0, 'Screensize'));
title('LF/HF'); xlabel('Permited error (\%)'); ylabel('Number of correct cases')



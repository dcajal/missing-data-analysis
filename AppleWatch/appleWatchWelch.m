%% Welch's method results with random distributed errors

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
wdw = hamming(60*fs);
noverlap = 30*fs;
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
            
            
%             figure;
%             plot(tReference,mtReference,'-'); hold on;
%             plot(taw,mtAwIpfm,'--');
%             plot(taw,mtAwIncidences,'-*');
%             plot(taw,mtAwIterative,'-^');
%             legend('Reference','IPFM','Incidences','Iterative')
%             axis tight; set(gcf, 'Position', get(0, 'Screensize'));
%             pause;
            
            
            % Remove baseline
            [bb, aa] = butter(4, 0.04*2/fs, 'high');
            auxMtReference= filtfilt(bb, aa, 1000./mtReference);
            auxMtAwIpfm = filtfilt(bb, aa, 1000./mtAwIpfm);
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
            
            if ~isempty(auxMtAwIpfm)
                psdAwIpfm = pwelch(auxMtAwIpfm,wdw,noverlap,f,fs);
                awIpfm{saveIndex}.hf = trapz(f(HFBand(1):HFBand(2)),psdAwIpfm(HFBand(1):HFBand(2)));
                awIpfm{saveIndex}.lf = trapz(f(LFBand(1):LFBand(2)),psdAwIpfm(LFBand(1):LFBand(2)));
                awIpfm{saveIndex}.lfn = awIpfm{saveIndex}.lf/(awIpfm{saveIndex}.lf+awIpfm{saveIndex}.hf);
                awIpfm{saveIndex}.lfhf = awIpfm{saveIndex}.lf/awIpfm{saveIndex}.hf;
            else
                awIpfm{saveIndex}.hf = nan;
                awIpfm{saveIndex}.lf = nan;
                awIpfm{saveIndex}.lfn = nan;
                awIpfm{saveIndex}.lfhf = nan;
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

            
%             figure;
%             plot(tReference,mtReference,'-'); hold on;
%             plot(taw,mtAwIpfm,'--');
%             plot(taw,mtAwIncidences,'-*');
%             plot(taw,mtAwIterative,'-^');
%             legend('Reference','IPFM','Incidences','Iterative')
%             axis tight; set(gcf, 'Position', get(0, 'Screensize'));
%             pause;
            
            
            % Remove baseline
            [bb, aa] = butter(4, 0.04*2/fs, 'high');
            auxMtReference= filtfilt(bb, aa, 1000./mtReference);
            auxMtAwIpfm = filtfilt(bb, aa, 1000./mtAwIpfm);
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

            if ~isempty(auxMtAwIpfm)
                psdAwIpfm = pwelch(auxMtAwIpfm,wdw,noverlap,f,fs);
                awIpfm{saveIndex}.hf = trapz(f(HFBand(1):HFBand(2)),psdAwIpfm(HFBand(1):HFBand(2)));
                awIpfm{saveIndex}.lf = trapz(f(LFBand(1):LFBand(2)),psdAwIpfm(LFBand(1):LFBand(2)));
                awIpfm{saveIndex}.lfn = awIpfm{saveIndex}.lf/(awIpfm{saveIndex}.lf+awIpfm{saveIndex}.hf);
                awIpfm{saveIndex}.lfhf = awIpfm{saveIndex}.lf/awIpfm{saveIndex}.hf;
            else
                awIpfm{saveIndex}.hf = nan;
                awIpfm{saveIndex}.lf = nan;
                awIpfm{saveIndex}.lfn = nan;
                awIpfm{saveIndex}.lfhf = nan;
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
errorThreshold = 0.001:0.001:0.3;

%% HF

ipfmCorrects = zeros(length(reference),numel(errorThreshold));
incidencesCorrects = zeros(length(reference),numel(errorThreshold));
iterativeCorrects = zeros(length(reference),numel(errorThreshold));
for kk=1:length(reference)
    ipfmCorrects(kk,:) = abs((awIpfm{kk}.hf - reference{kk}.hf)/reference{kk}.hf)<errorThreshold;
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
    ipfmCorrects(kk,:) = abs((awIpfm{kk}.lf - reference{kk}.lf)/reference{kk}.lf)<errorThreshold;
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
plot(100*errorThreshold,100*sum(ipfmCorrects,1)/numel(errorThreshold),'-');
plot(100*errorThreshold,100*sum(incidencesCorrects,1)/numel(errorThreshold),'-*');
plot(100*errorThreshold,100*sum(iterativeCorrects,1)/numel(errorThreshold),'-^');
legend('IPFM','Incidences','Iterative','Location','best')
axis tight; set(gcf, 'Position', get(0, 'Screensize'));
title('LF'); xlabel('Permited error (\%)'); ylabel('Correct cases (\%)')


%% LFn

ipfmCorrects = zeros(length(reference),numel(errorThreshold));
incidencesCorrects = zeros(length(reference),numel(errorThreshold));
iterativeCorrects = zeros(length(reference),numel(errorThreshold));
for kk=1:length(reference)
    ipfmCorrects(kk,:) = abs((awIpfm{kk}.lfn - reference{kk}.lfn)/reference{kk}.lfn)<errorThreshold;
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
    ipfmCorrects(kk,:) = abs((awIpfm{kk}.lfhf - reference{kk}.lfhf)/reference{kk}.lfhf)<errorThreshold;
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



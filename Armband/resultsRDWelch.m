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

debug = false;
segmentDurationSec = 3*60;
stepDurationSec = 1*60;
fs = 4;
nfft = 2^10;
wdw = hamming(60*fs);
noverlap = 30*fs;
f = linspace(0,0.4,nfft);
errorThreshold = 0.05;

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
            
            tholter = holterSegment(1):0.25:holterSegment(end);
            tarmband = armbandSegment(1):0.25:armbandSegment(end);
            
            
            % Reference (incidences)
            [tmholter,idsholter] = incidences(holterSegment,1);
            mtholter = ipfm(tmholter,idsholter,tholter);
            
            % Armband
            [tmarmband,idsarmband,~,tnarmband] = incidences(armbandSegment,1);
            mtarmbandipfm = ipfm(tmarmband,idsarmband,tarmband);
            mtarmbandincidences = ipfm(tnarmband,[],tarmband);
            mtarmbanditerative = ipfm(gapcorrector(armbandSegment),[],tarmband);
            
            if debug
                figure;
                plot(tholter,mtholter,'Linewidth',2); hold on
                plot(tarmband,mtarmbandipfm)
                plot(tarmband,mtarmbandincidences)
                plot(tarmband,mtarmbanditerative)
                legend('holter','ipfm','incidences','iterative')
            end

            % Remove baseline
            [bb, aa] = butter(4, 0.04*2/fs, 'high');
            auxmtholter= filtfilt(bb, aa, 1000./mtholter);
            auxmtarmbandipfm = filtfilt(bb, aa, 1000./mtarmbandipfm);
            auxmtarmbandincidences = filtfilt(bb, aa, 1000./mtarmbandincidences);
            auxmtarmbanditerative = filtfilt(bb, aa, 1000./mtarmbanditerative);
     
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % HF
            HFBand = [find(f>=0.15,1) find(f>=0.4,1)];
            LFBand = [find(f>=0.04,1) find(f>=0.15,1)];
            
            if ~isempty(auxmtholter)
                psdholter = pwelch(auxmtholter,wdw,noverlap,f,fs);
                holter{saveIndex}.hf = trapz(f(HFBand(1):HFBand(2)),psdholter(HFBand(1):HFBand(2)));
                holter{saveIndex}.lf = trapz(f(LFBand(1):LFBand(2)),psdholter(LFBand(1):LFBand(2)));
                holter{saveIndex}.lfn = holter{saveIndex}.lf/(holter{saveIndex}.lf+holter{saveIndex}.hf);
                holter{saveIndex}.lfhf = holter{saveIndex}.lf/holter{saveIndex}.hf;
            else
                holter{saveIndex}.hf = nan;
                holter{saveIndex}.lf = nan;
                holter{saveIndex}.lfn = nan;
                holter{saveIndex}.lfhf = nan;
            end
            
            if ~isempty(auxmtarmbandipfm)
                psdarmbandipfm = pwelch(auxmtarmbandipfm,wdw,noverlap,f,fs);
                armbandipfm{saveIndex}.hf = trapz(f(HFBand(1):HFBand(2)),psdarmbandipfm(HFBand(1):HFBand(2)));
                armbandipfm{saveIndex}.lf = trapz(f(LFBand(1):LFBand(2)),psdarmbandipfm(LFBand(1):LFBand(2)));
                armbandipfm{saveIndex}.lfn = armbandipfm{saveIndex}.lf/(armbandipfm{saveIndex}.lf+armbandipfm{saveIndex}.hf);
                armbandipfm{saveIndex}.lfhf = armbandipfm{saveIndex}.lf/armbandipfm{saveIndex}.hf;
            else
                armbandipfm{saveIndex}.hf = nan;
                armbandipfm{saveIndex}.lf = nan;
                armbandipfm{saveIndex}.lfn = nan;
                armbandipfm{saveIndex}.lfhf = nan;
            end
            
            if ~isempty(auxmtarmbandincidences)
                psdarmbandincidences = pwelch(auxmtarmbandincidences,wdw,noverlap,f,fs);
                armbandincidences{saveIndex}.hf = trapz(f(HFBand(1):HFBand(2)),psdarmbandincidences(HFBand(1):HFBand(2)));
                armbandincidences{saveIndex}.lf = trapz(f(LFBand(1):LFBand(2)),psdarmbandincidences(LFBand(1):LFBand(2)));
                armbandincidences{saveIndex}.lfn = armbandincidences{saveIndex}.lf/(armbandincidences{saveIndex}.lf+armbandincidences{saveIndex}.hf);
                armbandincidences{saveIndex}.lfhf = armbandincidences{saveIndex}.lf/armbandincidences{saveIndex}.hf;
            else
                armbandincidences{saveIndex}.hf = nan;
                armbandincidences{saveIndex}.lf = nan;
                armbandincidences{saveIndex}.lfn = nan;
                armbandincidences{saveIndex}.lfhf = nan;
            end
            
            if ~isempty(auxmtarmbanditerative)         
                psdarmbanditerative = pwelch(auxmtarmbanditerative,wdw,noverlap,f,fs);
                armbanditerative{saveIndex}.hf = trapz(f(HFBand(1):HFBand(2)),psdarmbanditerative(HFBand(1):HFBand(2)));
                armbanditerative{saveIndex}.lf = trapz(f(LFBand(1):LFBand(2)),psdarmbanditerative(LFBand(1):LFBand(2)));
                armbanditerative{saveIndex}.lfn = armbanditerative{saveIndex}.lf/(armbanditerative{saveIndex}.lf+armbanditerative{saveIndex}.hf);
                armbanditerative{saveIndex}.lfhf = armbanditerative{saveIndex}.lf/armbanditerative{saveIndex}.hf;
            else
                armbanditerative{saveIndex}.hf = nan;
                armbanditerative{saveIndex}.lf = nan;
                armbanditerative{saveIndex}.lfn = nan;
                armbanditerative{saveIndex}.lfhf = nan;
            end

            if debug
                figure;
                plot(f,psdholter,'Linewidth',2); hold on
                plot(f,psdarmbandipfm)
                plot(f,psdarmbandincidences)
                plot(f,psdarmbanditerative)
                legend('holter','ipfm','incidences','iterative')
            end
            
            saveIndex = saveIndex+1;
        end
    end             
end

%% HF

ipfmCorrects = zeros(length(holter),1);
incidencesCorrects = zeros(length(holter),1);
iterativeCorrects = zeros(length(holter),1);
for kk=1:length(holter)
    ipfmCorrects(kk) = abs((armbandipfm{kk}.hf - holter{kk}.hf)/holter{kk}.hf)<errorThreshold;
    incidencesCorrects(kk) = abs((armbandincidences{kk}.hf - holter{kk}.hf)/holter{kk}.hf)<errorThreshold;
    iterativeCorrects(kk) = abs((armbanditerative{kk}.hf - holter{kk}.hf)/holter{kk}.hf)<errorThreshold;
end
figure; hold on
p(1) = bar(1:length(holter),ipfmCorrects*3,'b');
patch([1 length(holter) length(holter) 1],[0 0 2 2],[1 1 1]);
p(2) = bar(1:length(holter),ipfmCorrects*2,'g');
patch([1 length(holter) length(holter) 1],[0 0 1 1],[1 1 1]);
p(3) = bar(1:length(holter),incidencesCorrects,'k');
legend(p,'IPFM','Incidences','Iterative','Location','bestoutside')
title('HF'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

fprintf('\n');
fprintf('----------------------------------------\n');
fprintf('----------------- HF -------------------\n');
fprintf('----------------------------------------\n');
fprintf('IPFM coverage: %i out of %i (%.2f%%)\n',sum(ipfmCorrects),numel(ipfmCorrects),sum(ipfmCorrects)/numel(ipfmCorrects)*100);
fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
fprintf('----------------------------------------\n');
fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
fprintf('Iterative coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(iterativeCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
fprintf('Incidences coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(incidencesCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
fprintf('----------------------------------------\n');

%% LF

ipfmCorrects = zeros(length(holter),1);
incidencesCorrects = zeros(length(holter),1);
iterativeCorrects = zeros(length(holter),1);
for kk=1:length(holter)
    ipfmCorrects(kk) = abs((armbandipfm{kk}.lf - holter{kk}.lf)/holter{kk}.lf)<errorThreshold;
    incidencesCorrects(kk) = abs((armbandincidences{kk}.lf - holter{kk}.lf)/holter{kk}.lf)<errorThreshold;
    iterativeCorrects(kk) = abs((armbanditerative{kk}.lf - holter{kk}.lf)/holter{kk}.lf)<errorThreshold;
end
figure; hold on
p(1) = bar(1:length(holter),ipfmCorrects*3,'b');
patch([1 length(holter) length(holter) 1],[0 0 2 2],[1 1 1]);
p(2) = bar(1:length(holter),ipfmCorrects*2,'g');
patch([1 length(holter) length(holter) 1],[0 0 1 1],[1 1 1]);
p(3) = bar(1:length(holter),incidencesCorrects,'k');
legend(p,'IPFM','Incidences','Iterative','Location','bestoutside')
title('LF'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

fprintf('\n');
fprintf('----------------------------------------\n');
fprintf('----------------- LF -------------------\n');
fprintf('----------------------------------------\n');
fprintf('IPFM coverage: %i out of %i (%.2f%%)\n',sum(ipfmCorrects),numel(ipfmCorrects),sum(ipfmCorrects)/numel(ipfmCorrects)*100);
fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
fprintf('----------------------------------------\n');
fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
fprintf('Iterative coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(iterativeCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
fprintf('Incidences coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(incidencesCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
fprintf('----------------------------------------\n');

%% LFn

ipfmCorrects = zeros(length(holter),1);
incidencesCorrects = zeros(length(holter),1);
iterativeCorrects = zeros(length(holter),1);
for kk=1:length(holter)
    ipfmCorrects(kk) = abs((armbandipfm{kk}.lfn - holter{kk}.lfn)/holter{kk}.lfn)<errorThreshold;
    incidencesCorrects(kk) = abs((armbandincidences{kk}.lfn - holter{kk}.lfn)/holter{kk}.lfn)<errorThreshold;
    iterativeCorrects(kk) = abs((armbanditerative{kk}.lfn - holter{kk}.lfn)/holter{kk}.lfn)<errorThreshold;
end
figure; hold on
p(1) = bar(1:length(holter),ipfmCorrects*3,'b');
patch([1 length(holter) length(holter) 1],[0 0 2 2],[1 1 1]);
p(2) = bar(1:length(holter),ipfmCorrects*2,'g');
patch([1 length(holter) length(holter) 1],[0 0 1 1],[1 1 1]);
p(3) = bar(1:length(holter),incidencesCorrects,'k');
legend(p,'IPFM','Incidences','Iterative','Location','bestoutside')
title('LFn'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

fprintf('\n');
fprintf('----------------------------------------\n');
fprintf('---------------- LFn -------------------\n');
fprintf('----------------------------------------\n');
fprintf('IPFM coverage: %i out of %i (%.2f%%)\n',sum(ipfmCorrects),numel(ipfmCorrects),sum(ipfmCorrects)/numel(ipfmCorrects)*100);
fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
fprintf('----------------------------------------\n');
fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
fprintf('Iterative coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(iterativeCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
fprintf('Incidences coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(incidencesCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
fprintf('----------------------------------------\n');

%%

ipfmCorrects = zeros(length(holter),1);
incidencesCorrects = zeros(length(holter),1);
iterativeCorrects = zeros(length(holter),1);
for kk=1:length(holter)
    ipfmCorrects(kk) = abs((armbandipfm{kk}.lfhf - holter{kk}.lfhf)/holter{kk}.lfhf)<errorThreshold;
    incidencesCorrects(kk) = abs((armbandincidences{kk}.lfhf - holter{kk}.lfhf)/holter{kk}.lfhf)<errorThreshold;
    iterativeCorrects(kk) = abs((armbanditerative{kk}.lfhf - holter{kk}.lfhf)/holter{kk}.lfhf)<errorThreshold;
end
figure; hold on
p(1) = bar(1:length(holter),ipfmCorrects*3,'b');
patch([1 length(holter) length(holter) 1],[0 0 2 2],[1 1 1]);
p(2) = bar(1:length(holter),ipfmCorrects*2,'g');
patch([1 length(holter) length(holter) 1],[0 0 1 1],[1 1 1]);
p(3) = bar(1:length(holter),incidencesCorrects,'k');
legend(p,'IPFM','Incidences','Iterative','Location','bestoutside')
title('LF/HF'); ylabel('Relative error'); axis tight; set(gcf,'position',[0,0,2000,1000]);

fprintf('\n');
fprintf('----------------------------------------\n');
fprintf('--------------- LF/HF ------------------\n');
fprintf('----------------------------------------\n');
fprintf('IPFM coverage: %i out of %i (%.2f%%)\n',sum(ipfmCorrects),numel(ipfmCorrects),sum(ipfmCorrects)/numel(ipfmCorrects)*100);
fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),numel(incidencesCorrects),sum(incidencesCorrects)/numel(incidencesCorrects)*100);
fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),numel(iterativeCorrects),sum(iterativeCorrects)/numel(iterativeCorrects)*100);
fprintf('----------------------------------------\n');
fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
fprintf('Iterative coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(iterativeCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
fprintf('Incidences coverage on IPFM failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~ipfmCorrects)),sum(~ipfmCorrects),sum(incidencesCorrects(~ipfmCorrects))/sum(~ipfmCorrects)*100);
fprintf('----------------------------------------\n');

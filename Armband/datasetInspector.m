clear
addpath('lib','database');
set(0,'defaultAxesFontName', 'Helvetica');
set(0,'defaultTextFontName', 'Helvetica');
set(0,'defaulttextInterpreter','latex')
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize', 20);
set(0,'defaultTextFontSize', 20);
set(0, 'defaultFigureRenderer', 'painters')

debug = false;
segmentDurationSec = 3*60;
stepDurationSec = 1*60;

% Load all files in database directory
dirlist = dir('database');
files = cell([1 length(dirlist)-2]);
for kk = 3:length(dirlist)
    files{kk-2} = dirlist(kk).name;
end

%%
kk = 229;
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

figure; plot(QRS2(2:end),diff(QRS2)*1000,'r'); hold on;

segmentBegin(artifactCounter>0) = [];
segmentEnd(artifactCounter>0) = [];
for jj = 1:length(segmentBegin)
    clean = QRS2(QRS2>segmentBegin(jj) & QRS2<segmentEnd(jj));
    plot(clean(2:end),diff(clean)*1000,'b');
end
function [ SupineECG, TiltECG ] = loadSimulationSubject( path, subject )
SIGNAL_DURATION = 120; % [s]
load(strcat(path,'/Simulation/database/',subject,'.mat'));

if SIGNAL_DURATION < 240
    tiltMarks(2) = tiltMarks(1) + SIGNAL_DURATION; %#ok<*NODEF>
    tiltMarks(4) = tiltMarks(3) + SIGNAL_DURATION;
end

% Split Supine-Tilt
% SupineECG.signal = Ecg.signal(tiltMarks(1)*Ecg.samplerate:tiltMarks(2)*Ecg.samplerate);
% SupineECG.t = 0:1/Ecg.samplerate:(length(SupineECG.signal)-1)/Ecg.samplerate;
SupineECG.tk = Ecg.qrs(find(Ecg.qrs>=tiltMarks(1),1):find(Ecg.qrs<=tiltMarks(2),1,'last'));
SupineECG.tk = SupineECG.tk-tiltMarks(1);
% SupineECG.samplerate = Ecg.samplerate;

% TiltECG.signal = Ecg.signal(tiltMarks(3)*Ecg.samplerate:tiltMarks(4)*Ecg.samplerate);
% TiltECG.t = 0:1/Ecg.samplerate:(length(TiltECG.signal)-1)/Ecg.samplerate;
TiltECG.tk = Ecg.qrs(find(Ecg.qrs>=tiltMarks(3),1):find(Ecg.qrs<=tiltMarks(4),1,'last'));
TiltECG.tk = TiltECG.tk-tiltMarks(3);
% TiltECG.samplerate = Ecg.samplerate;

SupineECG.tk = [2*SupineECG.tk(1)-SupineECG.tk(2); SupineECG.tk; 2*SupineECG.tk(end)-SupineECG.tk(end-1)];
SupineECG.tk = [2*SupineECG.tk(1)-SupineECG.tk(2); SupineECG.tk; 2*SupineECG.tk(end)-SupineECG.tk(end-1)];
TiltECG.tk = [2*TiltECG.tk(1)-TiltECG.tk(2); TiltECG.tk; 2*TiltECG.tk(end)-TiltECG.tk(end-1)];
TiltECG.tk = [2*TiltECG.tk(1)-TiltECG.tk(2); TiltECG.tk; 2*TiltECG.tk(end)-TiltECG.tk(end-1)];

end


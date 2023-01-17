%% Welch's method results with burst of missed beats

clear
addpath('lib','Simulation/database');
figurePresets
rng('default')

% Spectral analysis
windowSeconds = 60;
overlapSeconds = 30;
nfft = 2^10;
fs = 4;

% subject = {'01M'};
subject = {'01M' '02V' '04V' '05V' '06M' '07M' '11M' '13V' '17V'}; % Respiración en banda de HF
% subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '10V' '11M'...
%     '12V' '13V' '15V' '16V' '17V'};
burstDuration = [0 5 10 15 20];

fillGaps = 'iterativeNonLinear'; % 'ipfm' 'iterative' 'iterativeNonLinear'
nRealizations = 10;
displayCounter = 1;
useIPFM = true;
splineOrder = 4;

plotflag = false;
resultsSupineIPFM = cell(length(subject)*nRealizations,length(burstDuration));
resultsTiltIPFM = cell(length(subject)*nRealizations,length(burstDuration));
for kk = 1:length(subject)
    for ll = 1:nRealizations
        for jj = 1:length(burstDuration)
            fprintf('Computing realization %i of Subject %s with %i burst duration (%i/%i)...',...
                ll,subject{kk},burstDuration(jj),displayCounter,...
                length(subject)*nRealizations*length(burstDuration));
            [SupineECG, TiltECG] = loadSimulationSubject(pwd,subject{kk});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Pulse deletion
            if burstDuration(jj) > 0
                SupineECG.tk(find(SupineECG.tk<=30+(60-burstDuration(jj))/(nRealizations-1)*(ll-1),1,'Last'):...
                    find(SupineECG.tk>=30+(60-burstDuration(jj))/(nRealizations-1)*(ll-1)+burstDuration(jj),1)) = [];
                TiltECG.tk(find(TiltECG.tk<=30+(60-burstDuration(jj))/(nRealizations-1)*(ll-1),1,'Last'):...
                    find(TiltECG.tk>=30+(60-burstDuration(jj))/(nRealizations-1)*(ll-1)+burstDuration(jj),1)) = [];

                switch fillGaps
                    case 'iterativeNonLinear'
                        SupineECG.tn = gapcorrectorNonLinear(SupineECG.tk); SupineECG.ids = [];
                        TiltECG.tn = gapcorrectorNonLinear(TiltECG.tk); TiltECG.ids = [];
                    case 'iterative'
                        SupineECG.tn = gapcorrector(SupineECG.tk); SupineECG.ids = [];
                        TiltECG.tn = gapcorrector(TiltECG.tk); TiltECG.ids = [];
                    case 'incidences'
                        [~,~,~,SupineECG.tn] = incidences(SupineECG.tk,1); SupineECG.ids = [];
                        [~,~,~,TiltECG.tn] = incidences(TiltECG.tk,1); TiltECG.ids = [];
                    case 'ipfm'
                        % tn is actually tk. Only ids are computed
                        [SupineECG.tn,SupineECG.ids] = incidences(SupineECG.tk,1);
                        [TiltECG.tn,TiltECG.ids] = incidences(TiltECG.tk,1);
                end
            else
                SupineECG.tn = SupineECG.tk; SupineECG.ids = [];
                TiltECG.tn = TiltECG.tk; TiltECG.ids = [];
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % IPFM
            SupineIPFM.t = 0:1/fs:120;
            if useIPFM
                SupineIPFM.sp = mti(SupineECG.tn,SupineECG.ids,splineOrder);
                SupineIPFM.ihr = spval(SupineIPFM.sp,SupineIPFM.t);
            else
                error('Using filling gaps without IPFM');
                SupineIPFM.ihr = interp1(SupineECG.tn(2:end),1./diff(SupineECG.tn),SupineIPFM.t,'spline');   
            end

            TiltIPFM.t = SupineIPFM.t;
            if useIPFM
                TiltIPFM.sp = mti(TiltECG.tn,TiltECG.ids,splineOrder);
                TiltIPFM.ihr = spval(TiltIPFM.sp,TiltIPFM.t);
            else
                TiltIPFM.ihr = interp1(TiltECG.tn(2:end),1./diff(TiltECG.tn),TiltIPFM.t,'spline');
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Welch
            SupineIPFM.ihr = SupineIPFM.ihr*1000; % [mHz]
            TiltIPFM.ihr = TiltIPFM.ihr*1000; % [mHz]
            
            wdw = hamming(windowSeconds*fs);
            noverlap = overlapSeconds*fs;

            SupineIPFM.f = linspace(0,0.4,nfft);
            TiltIPFM.f = linspace(0,0.4,nfft);

            [bb, aa] = butter(4, 0.04*2/fs, 'high');
            auxSupine = filtfilt(bb, aa, SupineIPFM.ihr);
            auxTilt = filtfilt(bb, aa, TiltIPFM.ihr);

            SupineIPFM.pxx = pwelch(auxSupine,wdw,noverlap,SupineIPFM.f,fs);
            TiltIPFM.pxx = pwelch(auxTilt,wdw,noverlap,TiltIPFM.f,fs);
            clear wdw noverlap bb aa auxTilt auxSupine
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plots

            if plotflag
                if burstDuration(jj)==0
                    figure('DefaultAxesFontSize',14,'units','normalized','outerposition',[0 0 1 1]); 
                    subplot(211)
                    plot(SupineIPFM.t,SupineIPFM.ihr,'LineWidth',2); hold on
                    xlabel('Time [s]')
                    ylabel('IHR [mHz]')
                    title(subject{kk})
                    subplot(212)
                    plot(SupineIPFM.f,SupineIPFM.pxx,'LineWidth',2); hold on
                    xlabel('Frequency [Hz]')
                    ylabel('PSD_{Welch} [mHz^2]','interpreter','tex')
                    axis tight
                else
                    subplot(211)
                    plot(SupineIPFM.t,SupineIPFM.ihr,'LineWidth',1); hold on
                    subplot(212)
                    plot(SupineIPFM.f,SupineIPFM.pxx,'LineWidth',1); hold on
                end
                if burstDuration(jj)==20
                    legend('0','5','10','15','20');
        %             set(gcf,'position',[0,0,1000,300])
        %             set(gcf, 'Position', get(0, 'Screensize'));
                    pause
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            resultsSupineIPFM{(ll-1)*length(subject)+kk,jj} = freqind(SupineIPFM.pxx,SupineIPFM.f);
            resultsTiltIPFM{(ll-1)*length(subject)+kk,jj} = freqind(TiltIPFM.pxx,TiltIPFM.f);
            
            fprintf('Done\n');  
            displayCounter = displayCounter+1;
        end
    end
end

clear SupineECG SupineIPFM TiltECG TiltIPFM jj kk fs nfft overlapSeconds windowSeconds

%% Degradation results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Degradation results: Frequency indexes (Welch), error bursts, using %s method\n',fillGaps);
disp('Measure          Burst duration (seconds)');
fprintf('                 '); fprintf('%i                   ',burstDuration(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

error = computeError(resultsSupineIPFM, resultsTiltIPFM, burstDuration, 'LF');
% ylabel('P_{LF} [mHz^2]','interpreter','tex')
fprintf('P_LF           '); fprintf('%.2f (%.2f-%.2f)  & ',error); fprintf('\n')


error = computeError(resultsSupineIPFM, resultsTiltIPFM, burstDuration, 'HF');
% ylabel('P_{HF} [mHz^2]','interpreter','tex')
fprintf('P_HF           '); fprintf('%.2f (%.2f-%.2f)  & ',error); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')

clear error
%% Sympathovagal balance results
 
fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('p-values (supine/tilt groups): Frequency indexes (Welch), error bursts, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i          ',burstDuration); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

% figure(1)
% significance = twoGroupsDegradation(resultsSupineIPFM, resultsTiltIPFM, burstDuration, 'LF');
% ylabel('P_{LF} [mHz^2]','interpreter','tex')
% fprintf('P_LF           '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

figure(2)
significance = twoGroupsDegradation(resultsSupineIPFM, resultsTiltIPFM, burstDuration, 'HF');
ylabel('P_{HF} [mHz^2]','interpreter','tex')
fprintf('P_HF           '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

figure(3)
significance = twoGroupsDegradation(resultsSupineIPFM, resultsTiltIPFM, burstDuration, 'LFn');
ylabel('P_{LFn}','interpreter','tex')
fprintf('P_LFn          '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

figure(4)
significance = twoGroupsDegradation(resultsSupineIPFM, resultsTiltIPFM, burstDuration, 'LFHF');
ylabel('P_{LF}/P_{HF}','interpreter','tex')
fprintf('P_LF/P_HF      '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

clear significance

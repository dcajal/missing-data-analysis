%% Lomb's method results with random distributed errors

clear
addpath('lib','Simulation/database');
figurePresets
rng('default')

% Spectral analysis
windowSeconds = 60;
overlapSeconds = 30;
nfft = 2^10;

% subject = {'01M'};
subject = {'01M' '02V' '04V' '05V' '06M' '07M' '11M' '13V' '17Vmean'}; % Respiración en banda de HF
% subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '10V' '11M'...
%     '12V' '13V' '15V' '16V' '17V'};
% deletionProbability = 0:0.05:0.25;
deletionProbability = [0 0.05 0.15 0.25 0.35];

fillGaps = 'iterativeNonLinear'; % 'removeOutliers' 'iterative' 'iterativeNonLinear'
nRealizations = 10;
displayCounter = 1;
plotflag = false;

resultsSupineLomb = cell(length(subject)*nRealizations,length(deletionProbability));
resultsTiltLomb = cell(length(subject)*nRealizations,length(deletionProbability));
for kk = 1:length(subject)
    for ll = 1:nRealizations
        for jj = 1:length(deletionProbability)
            fprintf('Computing realization %i of Subject %s with %i%% deletion probability (%i/%i)...',...
                ll,subject{kk},deletionProbability(jj)*100,displayCounter,...
                length(subject)*nRealizations*length(deletionProbability));
            [SupineECG, TiltECG] = loadSimulationSubject(pwd,subject{kk});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Pulse deletion
            if deletionProbability(jj) > 0
                supineIndexes = binornd(1,deletionProbability(jj),length(SupineECG.tk),1)>0.5;
                supineIndexes(1:3) = 0;
                supineIndexes(end-3:end) = 0;
                tiltIndexes = binornd(1,deletionProbability(jj),length(TiltECG.tk),1)>0.5;
                tiltIndexes(1:3) = 0;
                tiltIndexes(end-3:end) = 0;

                SupineECG.tk(supineIndexes) = [];
                TiltECG.tk(tiltIndexes) = [];
                clear supineIndexes tiltIndexes

                switch fillGaps
                    case 'iterativeNonLinear'
                        SupineECG.tn = gapcorrectorNonLinear(SupineECG.tk); SupineECG.ids = [];
                        TiltECG.tn = gapcorrectorNonLinear(TiltECG.tk); TiltECG.ids = [];
                    case 'iterative'
                        SupineECG.tn = gapcorrector(SupineECG.tk);
                        TiltECG.tn = gapcorrector(TiltECG.tk);
                    case 'incidences'
                        [~,~,~,SupineECG.tn] = incidences(SupineECG.tk);
                        [~,~,~,TiltECG.tn] = incidences(TiltECG.tk);
                    case 'removeOutliers'
                        SupineECG.tn = SupineECG.tk;
                        TiltECG.tn = TiltECG.tk;
                end
            else
                SupineECG.tn = SupineECG.tk;
                TiltECG.tn = TiltECG.tk;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Lomb periodogram
            % Supine
            SupineLomb.f = linspace(0,0.4,nfft);
            SupineLomb.dtk = diff(SupineECG.tn);
            SupineLomb.tk = SupineECG.tn(2:end);

            % Delete peaks
            if strcmp(fillGaps,'removeOutliers')
                threshold = computeThreshold(SupineLomb.dtk);
                SupineLomb.tk(SupineLomb.dtk>threshold) = [];
                SupineLomb.dtk(SupineLomb.dtk>threshold) = [];
            end
            
            printAux = SupineLomb.dtk; 
            SupineLomb.hr = 1000./SupineLomb.dtk;                         
            SupineLomb.hr = detrend(SupineLomb.hr);
            clear threshold

            SupineLomb.pxx = plombav(SupineLomb.tk,SupineLomb.hr,SupineLomb.f,windowSeconds,overlapSeconds);

            [bb, aa] = butter(15, 0.04*2, 'high');
            h = freqz(bb,aa,nfft);
            SupineLomb.pxx = SupineLomb.pxx.*abs(h);
            SupineLomb.pxx = SupineLomb.pxx;
            clear h bb aa

            % Tilt
            TiltLomb.f = linspace(0,0.4,nfft);
            TiltLomb.dtk = diff(TiltECG.tn);
            TiltLomb.tk = TiltECG.tn(2:end);

            % Delete peaks
            if strcmp(fillGaps,'removeOutliers')
                threshold = computeThreshold(TiltLomb.dtk);
                TiltLomb.tk(TiltLomb.dtk>threshold) = [];
                TiltLomb.dtk(TiltLomb.dtk>threshold) = [];
            end
            TiltLomb.hr = 1000./TiltLomb.dtk;
            TiltLomb.hr = detrend(TiltLomb.hr);
            clear threshold

            TiltLomb.pxx = plombav(TiltLomb.tk,TiltLomb.hr,TiltLomb.f,windowSeconds,overlapSeconds);

            [bb, aa] = butter(15, 0.04*2, 'high');
            h = freqz(bb,aa,nfft);
            TiltLomb.pxx = TiltLomb.pxx.*abs(h);
            TiltLomb.pxx = TiltLomb.pxx;
            clear h bb aa        

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plots

            if plotflag
                if deletionProbability(jj)==0
                    figure('DefaultAxesFontSize',14,'units','normalized','outerposition',[0 0 1 1]); 
                    subplot(211)
                    stem(SupineLomb.tk,printAux,'LineWidth',2); hold on
                    ylabel('NN [s]');
                    xlabel('Time [s]');
                    title(subject{kk})
                    subplot(212)
                    plot(SupineLomb.f,SupineLomb.pxx,'LineWidth',2); hold on
                    xlabel('Frequency [Hz]')
                    ylabel('PSD_{Lomb} [Hz^2]','interpreter','tex')
                    axis tight
                else
                    subplot(211)
                    stem(SupineLomb.tk,printAux,'LineWidth',1); hold on
                    subplot(212)
                    plot(SupineLomb.f,SupineLomb.pxx,'LineWidth',1);
                end
                if deletionProbability(jj)==0.35
                    legend('0%','5%','15%','25%','35%');
        %             set(gcf,'position',[0,0,1000,300])
        %             set(gcf, 'Position', get(0, 'Screensize'));
                    pause
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            resultsSupineLomb{(ll-1)*length(subject)+kk,jj} = freqind(SupineLomb.pxx,SupineLomb.f);
            resultsTiltLomb{(ll-1)*length(subject)+kk,jj} = freqind(TiltLomb.pxx,TiltLomb.f);
            
            fprintf('Done\n');  
            displayCounter = displayCounter+1;
        end
    end
end

clear SupineECG SupineLomb TiltECG TiltLomb jj kk fs nfft overlapSeconds windowSeconds

%% Reference values

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('           Reference values\n');
disp('Measure');
fprintf('---------------------------------------------------------------------------------------\n')

reference = getReferenceValues(resultsSupineLomb, resultsTiltLomb, 'LF');
fprintf('P_LF             '); fprintf('%.2f (%.2f-%.2f)       ',reference); fprintf('\n')

reference = getReferenceValues(resultsSupineLomb, resultsTiltLomb, 'HF');
fprintf('P_HF            '); fprintf('%.2f (%.2f-%.2f)       ',reference); fprintf('\n')


fprintf('---------------------------------------------------------------------------------------\n')


%% Degradation results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Degradation results: Frequency indexes (Lomb), random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                '); fprintf('%i                      ',100*deletionProbability(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')


error = computeError(resultsSupineLomb, resultsTiltLomb, 100*deletionProbability, 'LF');
% ylabel('P_{LF} [mHz^2]','interpreter','tex')
fprintf('P_LF            '); fprintf('%.2f (%.2f -- %.2f) &  ',error); fprintf('\n')


error = computeError(resultsSupineLomb, resultsTiltLomb, 100*deletionProbability, 'HF');
% ylabel('P_{HF} [mHz^2]','interpreter','tex')
fprintf('P_HF            '); fprintf('%.2f (%.2f -- %.2f) &  ',error); fprintf('\n')


fprintf('---------------------------------------------------------------------------------------\n')

clear error

%% Sympathovagal balance results
 
fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('p-values (supine/tilt groups): Frequency indexes (Lomb), random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                '); fprintf('%i          ',deletionProbability(2:end)*100); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

% figure(1)
% significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, 100*deletionProbability, 'LF');
% ylabel('P_{LF} [mHz^2]','interpreter','tex')
% fprintf('P_LF           '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

figure(2)
significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, 100*deletionProbability, 'HF');
ylabel('P_{HF} [mHz^2]','interpreter','tex')
fprintf('P_HF           '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

figure(3)
significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, 100*deletionProbability, 'LFn');
ylabel('P_{LFn}','interpreter','tex')
fprintf('P_LFn          '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

figure(4)
significance = twoGroupsDegradation(resultsSupineLomb, resultsTiltLomb, 100*deletionProbability, 'LFHF');
ylabel('P_{LF}/P_{HF}','interpreter','tex')
fprintf('P_LF/P_HF      '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

clear significance
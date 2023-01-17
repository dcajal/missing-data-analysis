%% Poincare results with burst of missed beats

clear
addpath('lib','Simulation/database');
figurePresets
rng('default')

subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '10V' '11M'...
    '12V' '13V' '15V' '16V' '17V'};
burstDuration = [0 5 10 15 20];

fillGaps = 'iterativeNonLinear'; % 'removeOutliers' 'iterative' 'iterativeNonLinear'
detectGaps = false;
nRealizations = 10;
displayCounter = 1;

resultsSupineLPP = cell(length(subject)*nRealizations,length(burstDuration));
resultsTiltLPP = cell(length(subject)*nRealizations,length(burstDuration));
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
                        SupineECG.tn = gapcorrectorNonLinear(SupineECG.tk);
                        TiltECG.tn = gapcorrectorNonLinear(TiltECG.tk);
                    case 'iterative'
                        SupineECG.tn = gapcorrector(SupineECG.tk);
                        TiltECG.tn = gapcorrector(TiltECG.tk);
                    case 'incidences'
                        [~,~,~,SupineECG.tn] = incidences(SupineECG.tk);
                        [~,~,~,TiltECG.tn] = incidences(TiltECG.tk);
                    case 'noPreproc'
                        SupineECG.tn = SupineECG.tk;
                        TiltECG.tn = TiltECG.tk;
                    case 'removeOutliers'
                        SupineECG.tn = SupineECG.tk;
                        TiltECG.tn = TiltECG.tk;
                        detectGaps = true;
                end
            else
                SupineECG.tn = SupineECG.tk;
                TiltECG.tn = TiltECG.tk;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Lagged Poincare Plots
            SupineLPP = lpp(SupineECG.tn,1,detectGaps);
            TiltLPP = lpp(TiltECG.tn,1,detectGaps);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

            % Save results
            resultsSupineLPP{(ll-1)*length(subject)+kk,jj} = SupineLPP;
            resultsTiltLPP{(ll-1)*length(subject)+kk,jj} = TiltLPP;
            
            fprintf('Done\n');  
            displayCounter = displayCounter+1;
        end
    end
end

%% Degradation results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Degradation results: LPP indexes, error bursts, using %s method\n',fillGaps);
disp('Measure          Burst duration (seconds)');
fprintf('                 '); fprintf('%i                      ',burstDuration(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')


error = computeError(resultsSupineLPP, resultsTiltLPP, burstDuration, 'SD1');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('SD1 [ms]','interpreter','tex')
fprintf('SD1              '); fprintf('%.2f (%.2f -- %.2f) &  ',error); fprintf('\n')

error = computeError(resultsSupineLPP, resultsTiltLPP, burstDuration, 'SD2');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('SD2 [ms]','interpreter','tex')
fprintf('SD2              '); fprintf('%.2f (%.2f -- %.2f) &  ',error); fprintf('\n')

error = computeError(resultsSupineLPP, resultsTiltLPP, burstDuration, 'Md');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('Md   [ms]','interpreter','tex')
fprintf('Md               '); fprintf('%.2f (%.2f -- %.2f) &  ',error); fprintf('\n')

error = computeError(resultsSupineLPP, resultsTiltLPP, burstDuration, 'Sd');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('Sd   [ms]','interpreter','tex')
fprintf('Sd               '); fprintf('%.2f (%.2f -- %.2f) &  ',error); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')


%% Sympathovagal balance results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('p-values (supine/tilt groups): Poincare metrics, random distributed errors, using %s method\n',fillGaps);
disp('Measure          Burst duration (seconds)');
fprintf('                           '); fprintf('%i          ',burstDuration(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')


figure(1)
significance = twoGroupsDegradation(resultsSupineLPP, resultsTiltLPP, burstDuration, 'SD1');
fprintf('SD1           '); fprintf('%.3f       ',significance); fprintf('\n')

figure(2)
significance = twoGroupsDegradation(resultsSupineLPP, resultsTiltLPP, burstDuration, 'SD2');
fprintf('SD2           '); fprintf('%.3f       ',significance); fprintf('\n')

figure(3)
significance = twoGroupsDegradation(resultsSupineLPP, resultsTiltLPP, burstDuration, 'SD12');
fprintf('SD12          '); fprintf('%.3f       ',significance); fprintf('\n')

figure(4)
significance = twoGroupsDegradation(resultsSupineLPP, resultsTiltLPP, burstDuration, 'S');
fprintf('S             '); fprintf('%.3f       ',significance); fprintf('\n')

figure(5)
significance = twoGroupsDegradation(resultsSupineLPP, resultsTiltLPP, burstDuration, 'Md');
fprintf('Md            '); fprintf('%.3f       ',significance); fprintf('\n')

figure(6)
significance = twoGroupsDegradation(resultsSupineLPP, resultsTiltLPP, burstDuration, 'Sd');
fprintf('Sd            '); fprintf('%.3f       ',significance); fprintf('\n')

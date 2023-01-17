%% Poincare results with random distributed errors

clear
addpath('lib','Simulation/database');
figurePresets
rng('default')

subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '10V' '11M'...
    '12V' '13V' '15V' '16V' '17V'};
% deletionProbability = 0:0.05:0.25;
deletionProbability = [0 0.05 0.15 0.25 0.35];

fillGaps = 'iterativeNonLinear'; % 'removeOutliers' 'iterative' 'iterativeNonLinear'
detectGaps = false; % Only with fillGaps = 'remove outliers'
nRealizations = 10;
displayCounter = 1;

resultsSupineLPP = cell(length(subject)*nRealizations,length(deletionProbability));
resultsTiltLPP = cell(length(subject)*nRealizations,length(deletionProbability));
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

%% Reference values

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('           Reference values\n');
disp('Measure');
fprintf('---------------------------------------------------------------------------------------\n')

reference = getReferenceValues(resultsSupineLPP, resultsTiltLPP, 'SD1');
fprintf('SD1             '); fprintf('%.2f (%.2f-%.2f)       ',reference); fprintf('\n')

reference = getReferenceValues(resultsSupineLPP, resultsTiltLPP, 'SD2');
fprintf('SD2            '); fprintf('%.2f (%.2f-%.2f)       ',reference); fprintf('\n')

reference = getReferenceValues(resultsSupineLPP, resultsTiltLPP, 'Md');
fprintf('Md            '); fprintf('%.2f (%.2f-%.2f)       ',reference); fprintf('\n')

reference = getReferenceValues(resultsSupineLPP, resultsTiltLPP, 'Sd');
fprintf('Sd            '); fprintf('%.2f (%.2f-%.2f)       ',reference); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')


%% Degradation results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Degradation results: LPP indexes, random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i                     ',100*deletionProbability(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

error = computeError(resultsSupineLPP, resultsTiltLPP, 100*deletionProbability, 'SD1');
% ylabel('SD1 [ms]','interpreter','tex')
fprintf('SD1             '); fprintf('%.2f (%.2f-%.2f) &  ',error); fprintf('\n')

error = computeError(resultsSupineLPP, resultsTiltLPP, 100*deletionProbability, 'SD2');
% ylabel('SD2 [ms]','interpreter','tex')
fprintf('SD2             '); fprintf('%.2f (%.2f-%.2f) &  ',error); fprintf('\n')
 
error = computeError(resultsSupineLPP, resultsTiltLPP, 100*deletionProbability, 'Md');
% ylabel('Md   [ms]','interpreter','tex')
fprintf('Md              '); fprintf('%.2f (%.2f-%.2f) &  ',error); fprintf('\n')

error = computeError(resultsSupineLPP, resultsTiltLPP, 100*deletionProbability, 'Sd');
% ylabel('Sd   [ms]','interpreter','tex')
fprintf('Sd              '); fprintf('%.2f (%.2f-%.2f) &  ',error); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')


%% Sympathovagal balance results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('p-values (supine/tilt groups): Poincare metrics, random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                          '); fprintf('%i          ',deletionProbability(2:end)*100); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')


figure(1)
significance = twoGroupsDegradation(resultsSupineLPP, resultsTiltLPP, 100*deletionProbability, 'SD1');
fprintf('SD1           '); fprintf('%.3f       ',significance); fprintf('\n')

figure(2)
significance = twoGroupsDegradation(resultsSupineLPP, resultsTiltLPP, 100*deletionProbability, 'SD2');
fprintf('SD2           '); fprintf('%.3f       ',significance); fprintf('\n')

figure(3)
significance = twoGroupsDegradation(resultsSupineLPP, resultsTiltLPP, 100*deletionProbability, 'SD12');
fprintf('SD12          '); fprintf('%.3f       ',significance); fprintf('\n')

figure(4)
significance = twoGroupsDegradation(resultsSupineLPP, resultsTiltLPP, 100*deletionProbability, 'S');
fprintf('S             '); fprintf('%.3f       ',significance); fprintf('\n')

figure(5)
significance = twoGroupsDegradation(resultsSupineLPP, resultsTiltLPP, 100*deletionProbability, 'Md');
fprintf('Md            '); fprintf('%.3f       ',significance); fprintf('\n')

figure(6)
significance = twoGroupsDegradation(resultsSupineLPP, resultsTiltLPP, 100*deletionProbability, 'Sd');
fprintf('Sd            '); fprintf('%.3f       ',significance); fprintf('\n')

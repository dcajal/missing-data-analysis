%% Time-domain results with random distributed errors

clear
addpath('lib','Simulation/database');
figurePresets

subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '11M'...
    '12V' '13V' '15V' '16V' '17V'};
deletionProbability = 0:0.05:0.25;
fillGaps = 'removeOutliers'; % 'noPreproc' 'removeOutliers' 'incidences' 'iterative' 'iterativeNonLinear'
detectGaps = false; % Only with fillGaps = 'none' (remove outliers)

resultsSupineTDP = cell(length(subject),length(deletionProbability));
resultsTiltTDP = cell(length(subject),length(deletionProbability));
for kk = 1:length(subject)
    for jj = 1:length(deletionProbability)
        [SupineECG, TiltECG] = loadSimulationSubject(pwd,subject{kk});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Pulse deletion
        if deletionProbability(jj) > 0
            rng('default')
    
            SupineECG.tk(binornd(1,deletionProbability(jj),length(SupineECG.tk),1)>0.5) = [];
            TiltECG.tk(binornd(1,deletionProbability(jj),length(TiltECG.tk),1)>0.5) = [];

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

        % Time domain parameters
        [SupineTDP] = tempind(SupineECG.tn,detectGaps);
        [TiltTDP] = tempind(TiltECG.tn,detectGaps);
        
        
        % Save results
        resultsSupineTDP{kk,jj} = SupineTDP;
        resultsTiltTDP{kk,jj} = TiltTDP;
             
    end
end

clear SupineECG SupineTDP TiltECG TiltTDP jj kk

%% Degradation results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Degradation results: Time indexes, random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i                     ',100*deletionProbability(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

error = computeError(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'MHR');
% ylabel('MHR [beats/min]','interpreter','tex')
fprintf('MHR             '); fprintf('%.2f (%.2f-%.2f)       ',error); fprintf('\n')

error = computeError(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'SDNN');
% ylabel('SDNN [ms]','interpreter','tex')
fprintf('SDNN            '); fprintf('%.2f (%.2f-%.2f)       ',error); fprintf('\n')

error = computeError(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'SDSD');
% ylabel('SDSD [ms]','interpreter','tex')
fprintf('SDSD            '); fprintf('%.2f (%.2f-%.2f)       ',error); fprintf('\n')

% error = computeError(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'RMSSD');
% % ylabel('RMSSD [ms]','interpreter','tex')
% fprintf('RMSSD            '); fprintf('%.2f (%.2f-%.2f)        ',error); fprintf('\n')

error = computeError(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'pNN50');
% ylabel('pNN50 [%]','interpreter','tex')
fprintf('pNN50           '); fprintf('%.2f (%.2f-%.2f)       ',error); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')


%% Sympathovagal balance results
 
fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('p-values (supine/tilt groups): Time indexes, random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i          ',deletionProbability(2:end)*100); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'MHR');
ylabel('MHR [beats/min]','interpreter','tex')
fprintf('MHR              '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'SDNN');
ylabel('SDNN [ms]','interpreter','tex')
fprintf('SDNN             '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'SDSD');
ylabel('SDSD [ms]','interpreter','tex')
fprintf('SDSD             '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

% significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'RMSSD');
% ylabel('RMSSD [ms]','interpreter','tex')
% fprintf('RMSSD            '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'pNN50');
ylabel('pNN50 [%]','interpreter','tex')
fprintf('pNN50            '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')

clear significance
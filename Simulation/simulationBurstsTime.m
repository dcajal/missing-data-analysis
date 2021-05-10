%% Time-domain results with burst of missed beats

clear
addpath('lib','Simulation/database');
figurePresets

subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '11M'...
    '12V' '13V' '15V' '16V' '17V'};
burstDuration = [0 10 20 30 40 50 60];
fillGaps = 'iterativeNonLinear'; % 'noPreproc' 'removeOutliers' 'incidences' 'iterative' 'iterativeNonLinear'
detectGaps = false; % Only with fillGaps = 'none' (remove outliers)

resultsSupineTDP = cell(length(subject),length(burstDuration));
resultsTiltTDP = cell(length(subject),length(burstDuration));
for kk = 1:length(subject)
    for jj = 1:length(burstDuration)
        [SupineECG, TiltECG] = loadSimulationSubject(pwd,subject{kk});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Pulse deletion
        if burstDuration(jj) > 0
            rng('default')
    
            SupineECG.tk(find(SupineECG.tk<=60-burstDuration(jj)/2,1,'Last'):...
                find(SupineECG.tk>=60+burstDuration(jj)/2,1)) = [];
            TiltECG.tk(find(TiltECG.tk<=60-burstDuration(jj)/2,1,'Last'):...
                find(TiltECG.tk>=60+burstDuration(jj)/2,1)) = [];

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
fprintf('Degradation results: Time indexes, error bursts, using %s method\n',fillGaps);
disp('Measure          Burst duration (seconds)');
fprintf('                 '); fprintf('%i                      ',burstDuration(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')


error = computeError(resultsSupineTDP, resultsTiltTDP, burstDuration, 'MHR');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('MHR [beats/min]','interpreter','tex')
fprintf('MHR              '); fprintf('%.2f (%.2f-%.2f)        ',error); fprintf('\n')

error = computeError(resultsSupineTDP, resultsTiltTDP, burstDuration, 'SDNN');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('SDNN [ms]','interpreter','tex')
fprintf('SDNN             '); fprintf('%.2f (%.2f-%.2f)        ',error); fprintf('\n')

error = computeError(resultsSupineTDP, resultsTiltTDP, burstDuration, 'SDSD');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('SDSD [ms]','interpreter','tex')
fprintf('SDSD             '); fprintf('%.2f (%.2f-%.2f)        ',error); fprintf('\n')

% error = computeError(resultsSupineTDP, resultsTiltTDP, burstDuration, 'RMSSD');
% % xlabel('Burst duration (s)','interpreter','tex')
% % ylabel('RMSSD [ms]','interpreter','tex')
% fprintf('RMSSD            '); fprintf('%.2f (%.2f-%.2f)        ',error); fprintf('\n')

error = computeError(resultsSupineTDP, resultsTiltTDP, burstDuration, 'pNN50');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('pNN50 [%]','interpreter','tex')
fprintf('pNN50            '); fprintf('%.2f (%.2f-%.2f)        ',error); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')


%% Sympathovagal balance results
 
fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('p-values (supine/tilt groups): Time indexes, error bursts, using %s method\n',fillGaps);
disp('Measure          Burst duration (seconds)');
fprintf('                 '); fprintf('%i          ',burstDuration(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, burstDuration, 'MHR');
ylabel('MHR [beats/min]','interpreter','tex')
fprintf('MHR              '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, burstDuration, 'SDNN');
ylabel('SDNN [ms]','interpreter','tex')
fprintf('SDNN             '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, burstDuration, 'SDSD');
ylabel('SDSD [ms]','interpreter','tex')
fprintf('SDSD             '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

% significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, burstDuration, 'RMSSD');
% ylabel('RMSSD [ms]','interpreter','tex')
% fprintf('RMSSD            '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, burstDuration, 'pNN50');
ylabel('pNN50 [%]','interpreter','tex')
fprintf('pNN50            '); fprintf('%.3f       ',significance(2:end)); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')

clear significance
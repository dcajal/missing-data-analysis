%% Time-domain results with burst of missed beats

clear
addpath('lib','Simulation/database');
figurePresets
rng('default')

subject = {'01M' '02V' '03V' '04V' '05V' '06M' '07M' '08M' '09M' '10V' '11M'...
    '12V' '13V' '15V' '16V' '17V'};
burstDuration = [0 5 10 15 20];

fillGaps = 'iterativeNonLinear'; % 'removeOutliers' 'iterative' 'iterativeNonLinear'
detectGaps = false; % Only with fillGaps = removeOutliers (automatic)
nRealizations = 10;
displayCounter = 1;

tic
resultsSupineTDP = cell(length(subject)*nRealizations,length(burstDuration));
resultsTiltTDP = cell(length(subject)*nRealizations,length(burstDuration));
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

            % Time domain parameters
            [SupineTDP] = tempind(SupineECG.tn,detectGaps);
            [TiltTDP] = tempind(TiltECG.tn,detectGaps);


            % Save results
            resultsSupineTDP{(ll-1)*length(subject)+kk,jj} = SupineTDP;
            resultsTiltTDP{(ll-1)*length(subject)+kk,jj} = TiltTDP;

            fprintf('Done\n');  
            displayCounter = displayCounter+1;
        end
    end
end
toc

clear SupineECG SupineTDP TiltECG TiltTDP jj kk

%% Degradation results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Degradation results: Time indexes, error bursts, using %s method\n',fillGaps);
disp('Measure       Burst duration (seconds)');
fprintf('              '); fprintf('%i                      ',burstDuration(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')


error = computeError(resultsSupineTDP, resultsTiltTDP, burstDuration, 'MHR');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('MHR [beats/min]','interpreter','tex')
fprintf('MHR           '); fprintf('%.2f (%.2f-%.2f)  & ',error); fprintf('\n')

error = computeError(resultsSupineTDP, resultsTiltTDP, burstDuration, 'SDNN');
% xlabel('Burst duration (s)','interpreter','tex')
% ylabel('SDNN [ms]','interpreter','tex')
fprintf('SDNN          '); fprintf('%.2f (%.2f-%.2f)  & ',error); fprintf('\n')

error = computeError(resultsSupineTDP, resultsTiltTDP, burstDuration, 'RMSSD');
% % xlabel('Burst duration (s)','interpreter','tex')
% % ylabel('RMSSD [ms]','interpreter','tex')
fprintf('RMSSD         '); fprintf('%.2f (%.2f-%.2f)  & ',error); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')


%% Sympathovagal balance results
 
fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('p-values (supine/tilt groups): Time indexes, error bursts, using %s method\n',fillGaps);
disp('Measure          Burst duration (seconds)');
fprintf('                 '); fprintf('%i          ',burstDuration(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

figure(1)
significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, burstDuration, 'MHR');
ylabel('MHR [beats/min]','interpreter','tex')
% printeps(1,'burst_ro_mhr');
fprintf('MHR              '); fprintf('%.3f       ',significance); fprintf('\n')

figure(2)
significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, burstDuration, 'SDNN');
ylabel('SDNN [ms]','interpreter','tex')
fprintf('SDNN             '); fprintf('%.3f       ',significance); fprintf('\n')

figure(3)
significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, burstDuration, 'RMSSD');
ylabel('RMSSD [ms]','interpreter','tex')
% printeps(3,'burst_ro_rmssd');
fprintf('RMSSD            '); fprintf('%.3f       ',significance); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')

clear significance
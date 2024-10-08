%% Time-domain results with random distributed errors

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

tic
resultsSupineTDP = cell(length(subject)*nRealizations,length(deletionProbability));
resultsTiltTDP = cell(length(subject)*nRealizations,length(deletionProbability));
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

%% Reference values

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('           Reference values\n');
disp('Measure');
fprintf('---------------------------------------------------------------------------------------\n')

reference = getReferenceValues(resultsSupineTDP, resultsTiltTDP, 'MHR');
fprintf('MHR             '); fprintf('%.2f (%.2f-%.2f)       ',reference); fprintf('\n')

reference = getReferenceValues(resultsSupineTDP, resultsTiltTDP, 'SDNN');
fprintf('SDNN            '); fprintf('%.2f (%.2f-%.2f)       ',reference); fprintf('\n')

reference = getReferenceValues(resultsSupineTDP, resultsTiltTDP, 'RMSSD');
fprintf('RMSSD            '); fprintf('%.2f (%.2f-%.2f)        ',reference); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')


%% Degradation results

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Degradation results: Time indexes, random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i                     ',100*deletionProbability(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

[error,MHRp,MHRrvalues] = computeError(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'MHR');
% ylabel('MHR [beats/min]','interpreter','tex')
fprintf('MHR             '); fprintf('%.2f (%.2f-%.2f)  & ',error); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f              & ',rvalues); fprintf('\n')


[error,SDNNp,SDNNrvalues] = computeError(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'SDNN');
% ylabel('SDNN [ms]','interpreter','tex')
fprintf('SDNN            '); fprintf('%.2f (%.2f-%.2f)  & ',error); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f              & ',rvalues); fprintf('\n')


[error,RMSSDp,RMSSDrvalues] = computeError(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'RMSSD');
% ylabel('RMSSD [ms]','interpreter','tex')
fprintf('RMSSD           '); fprintf('%.2f (%.2f-%.2f)  & ',error); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f              & ',rvalues); fprintf('\n')


fprintf('---------------------------------------------------------------------------------------\n')


%% Sympathovagal balance results
 
fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('p-values (supine/tilt groups): Time indexes, random distributed errors, using %s method\n',fillGaps);
disp('Measure          Deletion probability (%)');
fprintf('                 '); fprintf('%i          ',deletionProbability(2:end)*100); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

figure(1)
significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'MHR');
ylabel('MHR [beats/min]')
% printeps(1,'rd_il_mhr');
fprintf('MHR              '); fprintf('%.3f       ',significance); fprintf('\n')

figure(2)
significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'SDNN');
ylabel('SDNN [ms]')
fprintf('SDNN             '); fprintf('%.3f       ',significance); fprintf('\n')

figure(3)
significance = twoGroupsDegradation(resultsSupineTDP, resultsTiltTDP, 100*deletionProbability, 'RMSSD');
ylabel('RMSSD [ms]');
% printeps(3,'rd_il_rmssd');
fprintf('RMSSD            '); fprintf('%.3f       ',significance); fprintf('\n')

fprintf('---------------------------------------------------------------------------------------\n')

clear significance
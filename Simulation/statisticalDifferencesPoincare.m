%% Poincare

clear
addpath('lib','Simulation/database');
deletionProbability = [0 0.05 0.15 0.25 0.35];
category = 'BurstsPoincare';

load(strcat('Simulation/results/',category,'/removeOutliers.mat'))
ro = [resultsSupineLPP; resultsTiltLPP]; clear resultsSupineLPP resultsTiltLPP

load(strcat('Simulation/results/',category,'/iterative.mat'))
l = [resultsSupineLPP; resultsTiltLPP]; clear resultsSupineLPP resultsTiltLPP

load(strcat('Simulation/results/',category,'/iterativeNonLinear.mat'))
nl = [resultsSupineLPP; resultsTiltLPP]; clear resultsSupineLPP resultsTiltLPP


%% pvalues
metric = 'Sd';

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
disp('Measure          Deletion probability (%)');
fprintf('                '); fprintf('%i         ',100*deletionProbability(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

figure(1)
significance = twoGroupsDegradation(ro, l, 100*deletionProbability, metric);
fprintf('pvalues ro-l      '); fprintf('%.2f       ',significance(2:end)); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')

figure(2)
significance = twoGroupsDegradation(l, nl, 100*deletionProbability, metric);
fprintf('pvalues l-nl      '); fprintf('%.2f       ',significance(2:end)); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')

figure(3)
significance = twoGroupsDegradation(ro, nl, 100*deletionProbability, metric);
fprintf('pvalues ro-nl     '); fprintf('%.2f       ',significance(2:end)); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')


fprintf('---------------------------------------------------------------------------------------\n')

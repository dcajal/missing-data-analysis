%% Welch

clear
addpath('lib','Simulation/database');
deletionProbability = [0 0.05 0.15 0.25 0.35];
category = 'BurstsWelch';

load(strcat('Simulation/results/',category,'/ipfm.mat'))
ipfm = [resultsSupineIPFM; resultsTiltIPFM]; clear resultsSupineIPFM resultsTiltIPFM

load(strcat('Simulation/results/',category,'/iterative.mat'))
l = [resultsSupineIPFM; resultsTiltIPFM]; clear resultsSupineIPFM resultsTiltIPFM

load(strcat('Simulation/results/',category,'/iterativeNonLinear.mat'))
nl = [resultsSupineIPFM; resultsTiltIPFM]; clear resultsSupineIPFM resultsTiltIPFM


%% pvalues
metric = 'HF';

fprintf('\n'); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')
disp('Measure          Deletion probability (%)');
fprintf('                '); fprintf('%i         ',100*deletionProbability(2:end)); fprintf('\n')
fprintf('---------------------------------------------------------------------------------------\n')

figure(1)
significance = twoGroupsDegradation(ipfm, l, 100*deletionProbability, metric);
fprintf('pvalues ipfm-l    '); fprintf('%.2f       ',significance(2:end)); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')

figure(2)
significance = twoGroupsDegradation(l, nl, 100*deletionProbability, metric);
fprintf('pvalues l-nl      '); fprintf('%.2f       ',significance(2:end)); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')

figure(3)
significance = twoGroupsDegradation(ipfm, nl, 100*deletionProbability, metric);
fprintf('pvalues ipfm-nl   '); fprintf('%.2f       ',significance(2:end)); fprintf('\n')
% fprintf('rvalues         '); fprintf('%.2f       ',rvalues(2:end)); fprintf('\n')


fprintf('---------------------------------------------------------------------------------------\n')

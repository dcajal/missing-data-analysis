function [ errorOutput, pOutput, rvalues ] = computeError( resultsSupine, resultsTilt, deletionProbability, index)

nsubjects = size(resultsSupine,1);

aux1 = [];
significance = 1;
zvalue = [];
for jj=1:length(deletionProbability)
    for kk=1:nsubjects
        aux1(kk,jj) = resultsSupine{kk,jj}.(index); %#ok<AGROW,*SAGROW>
    end
    for kk=nsubjects+1:2*nsubjects
        aux1(kk,jj) = resultsTilt{kk-nsubjects,jj}.(index); %#ok<AGROW,*SAGROW>
    end
    if jj>1
        [significance(jj),~,zvalueAux] = signrank(aux1(:,1),aux1(:,jj),'method','approximate'); %#ok<AGROW>
        zvalue(jj) = zvalueAux.zval; %#ok<AGROW>
    end
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;

% Boxplot
% figure
% fancyBoxplot(aux2,aux1,deletionProbability,significance)

rvalues = abs(zvalue)/sqrt(size(aux1,1));

% error = abs(aux2-aux1); % Absolute
error = abs(aux2-aux1)./aux2*100; % Relative
% error = (aux2-aux1)./aux2*100; % Relative (with sign)
errorMedian = prctile(error,50);
firstQuartile = prctile(error,25);
thirdQuartile = prctile(error,75);

[p,h]=signrank(error(:,2),zeros(size(error(:,2))),'method','approximate'); pOutput = p;
[p,h]=signrank(error(:,3),zeros(size(error(:,2))),'method','approximate'); pOutput = [pOutput p];
[p,h]=signrank(error(:,4),zeros(size(error(:,2))),'method','approximate'); pOutput = [pOutput p];
[p,h]=signrank(error(:,5),zeros(size(error(:,2))),'method','approximate'); pOutput = [pOutput p];

errorOutput = [];
for kk = 2:length(deletionProbability)
   errorOutput = [errorOutput errorMedian(kk) firstQuartile(kk) thirdQuartile(kk)];
end

end


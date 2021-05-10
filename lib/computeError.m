function [ output ] = computeError( resultsSupine, resultsTilt, deletionProbability, index)

nsubjects = size(resultsSupine,1);

aux1 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:nsubjects
        aux1(kk,jj) = resultsSupine{kk,jj}.(index); %#ok<AGROW,*SAGROW>
    end
    for kk=nsubjects+1:2*nsubjects
        aux1(kk,jj) = resultsTilt{kk-nsubjects,jj}.(index); %#ok<AGROW,*SAGROW>
    end
    significance(jj) = signrank(aux1(:,1),aux1(:,jj)); %#ok<AGROW>
end
aux2 = repmat(aux1(:,1),1,length(deletionProbability));
aux1(:,1) = nan;
% fancyBoxplot(aux2,aux1,deletionProbability,significance)
% ylabel(index,'interpreter','tex');

error = abs(aux2-aux1);
errorMedian = prctile(error,50);
firstQuartile = prctile(error,25);
thirdQuartile = prctile(error,75);

output = [];
for kk = 2:length(deletionProbability)
   output = [output errorMedian(kk) firstQuartile(kk) thirdQuartile(kk)];
end

end


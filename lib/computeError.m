function [ RMSE ] = computeError( resultsSupine, resultsTilt, deletionProbability, index, normalize)
if nargin < 5, normalize = false; end

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
fancyBoxplot(aux2,aux1,deletionProbability,significance)
ylabel(index,'interpreter','tex');

if normalize
    RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*nsubjects))./nanmean(aux2(:,1)); % nRMSE
else
    RMSE = sqrt(nansum((aux2-aux1).^2,1)./(2*nsubjects)); % RMSE
end

end


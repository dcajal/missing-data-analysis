function [ significance ] = twoGroupsDegradation( resultsSupine, resultsTilt, deletionProbability, index)
nsubjects = size(resultsSupine,1);

aux1 = [];
aux2 = [];
significance = [];
for jj=1:length(deletionProbability)
    for kk=1:nsubjects
        aux1(kk,jj) = resultsSupine{kk,jj}.(index); %#ok<AGROW>
        aux2(kk,jj) = resultsTilt{kk,jj}.(index); %#ok<AGROW>
    end
    significance(jj) = signrank(aux1(:,jj),aux2(:,jj)); %#ok<AGROW>
end
% fancyBoxplot(aux1,aux2,deletionProbability,significance,true)
% ylabel(index,'interpreter','tex');

end


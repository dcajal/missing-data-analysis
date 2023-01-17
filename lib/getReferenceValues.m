function [ reference ] = getReferenceValues( resultsSupine, resultsTilt, index )

nsubjects = size(resultsSupine,1);

aux1 = [];
for kk=1:nsubjects
    aux1(kk) = resultsSupine{kk,1}.(index); %#ok<AGROW,*SAGROW>
end
for kk=nsubjects+1:2*nsubjects
    aux1(kk) = resultsTilt{kk-nsubjects,1}.(index); %#ok<AGROW,*SAGROW>
end

referenceMedian = prctile(aux1,50);
firstQuartile = prctile(aux1,25);
thirdQuartile = prctile(aux1,75);

reference = [referenceMedian firstQuartile thirdQuartile];

end


function [] = groupDiscriminationAW(results, index)

nsubjects = size(results,2);

for kk=1:nsubjects
    reference(kk) = results{1,kk}.(index); %#ok<AGROW>
    OR(kk) = results{2,kk}.(index); %#ok<AGROW>
    L(kk) = results{3,kk}.(index); %#ok<AGROW>
    NL(kk) = results{4,kk}.(index); %#ok<AGROW>
end

significance(1) = signrank(reference(1:2:end), reference(2:2:end));
significance(2) = signrank(OR(1:2:end), OR(2:2:end));
significance(3) = signrank(L(1:2:end), L(2:2:end));
significance(4) = signrank(NL(1:2:end), NL(2:2:end));

fancyBoxplotAW(reference,OR,L,NL,significance); set(gcf,'position',[100,100,600,450]); box off

end


function [ output, p ] = computeWelchCoverage( input, errorThreshold, index )
reference = input(1,:);
M = input(2,:);
L = input(3,:);
NL = input(4,:);

nCases = length(reference);

MCorrects = zeros(nCases,numel(errorThreshold));
LCorrects = zeros(nCases,numel(errorThreshold));
NLCorrects = zeros(nCases,numel(errorThreshold));
MError = zeros(nCases,1);
LError = zeros(nCases,1);
NLError = zeros(nCases,1);
for kk=1:nCases
    MCorrects(kk,:) = abs((M{kk}.(index) - reference{kk}.(index))/reference{kk}.(index))<errorThreshold;
    LCorrects(kk,:) = abs((L{kk}.(index) - reference{kk}.(index))/reference{kk}.(index))<errorThreshold;
    NLCorrects(kk,:) = abs((NL{kk}.(index) - reference{kk}.(index))/reference{kk}.(index))<errorThreshold;
    
    % Absolute
%     MError(kk) = abs(M{kk}.(index) - reference{kk}.(index));
%     LError(kk) = abs(L{kk}.(index) - reference{kk}.(index));
%     NLError(kk) = abs(NL{kk}.(index) - reference{kk}.(index));
    
    % Relative
    MError(kk) = abs(M{kk}.(index) - reference{kk}.(index))/reference{kk}.(index);
    LError(kk) = abs(L{kk}.(index) - reference{kk}.(index))/reference{kk}.(index);
    NLError(kk) = abs(NL{kk}.(index) - reference{kk}.(index))/reference{kk}.(index);
end

if numel(errorThreshold) > 1
%     figure; 
    hold on
    p(1) = plot(errorThreshold*100,100*sum(MCorrects,1)/nCases,'b','Linewidth',1.5);
    p(2) = plot(errorThreshold*100,100*sum(LCorrects,1)/nCases,'g','Linewidth',1.5);
    p(3) = plot(errorThreshold*100,100*sum(NLCorrects,1)/nCases,'r','Linewidth',1.5);
%     title(index); 
    xlabel('Permitted error (%)'); ylabel('Correct cases (%)')
else
    fprintf('\n');
    fprintf('----------------------------------------\n');
    fprintf('---------------- %s -------------------\n',index);
    fprintf('----------------------------------------\n');
    fprintf('M coverage: %i out of %i (%.2f%%)\n',sum(MCorrects),nCases,sum(MCorrects)/nCases*100);
    fprintf('L coverage: %i out of %i (%.2f%%)\n',sum(LCorrects),nCases,sum(LCorrects)/nCases*100);
    fprintf('NL coverage: %i out of %i (%.2f%%)\n',sum(NLCorrects),nCases,sum(NLCorrects)/nCases*100);
    fprintf('----------------------------------------\n');
    fprintf('L coverage on M failures: %i out of %i (%.2f%%)\n',sum(LCorrects(~MCorrects)),sum(~MCorrects),sum(LCorrects(~MCorrects))/sum(~MCorrects)*100);
    fprintf('----------------------------------------\n');

    figure; hold on
    patch([0 nCases+1 nCases+1 0],[0 0 3 3],[1 1 1]);
    p(1) = bar(1:nCases,MCorrects*3,'g');
    patch([0 nCases+1 nCases+1 0],[0 0 2 2],[1 1 1]);
    p(2) = bar(1:nCases,LCorrects*2,'r');
    patch([0 nCases+1 nCases+1 0],[0 0 1 1],[1 1 1]);
    p(3) = bar(1:nCases,NLCorrects,'k');
    yticks([]);
    title(index); xlabel('Case'); ylabel('Is correct');
end
    axis tight;
    legend(p,'M','L','NL','Location','best');
    
    % Error
    output = [];
    
    errorMedian = prctile(MError,50);
    firstQuartile = prctile(MError,25);
    thirdQuartile = prctile(MError,75);
    output = [output errorMedian firstQuartile thirdQuartile];

    errorMedian = prctile(LError,50);
    firstQuartile = prctile(LError,25);
    thirdQuartile = prctile(LError,75);
    output = [output errorMedian firstQuartile thirdQuartile];

    errorMedian = prctile(NLError,50);
    firstQuartile = prctile(NLError,25);
    thirdQuartile = prctile(NLError,75);
    output = [output errorMedian firstQuartile thirdQuartile];

end


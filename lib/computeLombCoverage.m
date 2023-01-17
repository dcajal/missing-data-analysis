function [ output, p ] = computeLombCoverage( input, errorThreshold, index )
reference = input(1,:);
OR = input(2,:);
L = input(3,:);
NL = input(4,:);

nCases = length(reference);

ORCorrects = zeros(nCases,numel(errorThreshold));
LCorrects = zeros(nCases,numel(errorThreshold));
NLCorrects = zeros(nCases,numel(errorThreshold));

ORError = zeros(nCases,1);
LError = zeros(nCases,1);
NLError = zeros(nCases,1);
for kk=1:nCases
    ORCorrects(kk,:) = abs((OR{kk}.(index) - reference{kk}.(index))/reference{kk}.(index))<errorThreshold;
    LCorrects(kk,:) = abs((L{kk}.(index) - reference{kk}.(index))/reference{kk}.(index))<errorThreshold;
    NLCorrects(kk,:) = abs((NL{kk}.(index) - reference{kk}.(index))/reference{kk}.(index))<errorThreshold;
    
    % Absolute
%     ORError(kk) = abs(OR{kk}.(index) - reference{kk}.(index));
%     LError(kk) = abs(L{kk}.(index) - reference{kk}.(index));
%     NLError(kk) = abs(NL{kk}.(index) - reference{kk}.(index));
    
    % Relative
%     ORError(kk) = abs(OR{kk}.(index) - reference{kk}.(index))/reference{kk}.(index);
%     LError(kk) = abs(L{kk}.(index) - reference{kk}.(index))/reference{kk}.(index);
%     NLError(kk) = abs(NL{kk}.(index) - reference{kk}.(index))/reference{kk}.(index);
    
    % Relative (with sign)
    ORError(kk) = (OR{kk}.(index) - reference{kk}.(index))/reference{kk}.(index);
    LError(kk) = (L{kk}.(index) - reference{kk}.(index))/reference{kk}.(index);
    NLError(kk) = (NL{kk}.(index) - reference{kk}.(index))/reference{kk}.(index);
end

if numel(errorThreshold) > 1
%     figure; 
    hold on
    p(1) = plot(errorThreshold*100,100*sum(ORCorrects,1)/nCases,'b','Linewidth',1.5);
    p(2) = plot(errorThreshold*100,100*sum(LCorrects,1)/nCases,'g','Linewidth',1.5);
    p(3) = plot(errorThreshold*100,100*sum(NLCorrects,1)/nCases,'r','Linewidth',1.5);
%     title(index);
    xlabel('Permitted error (%)','interpreter','tex'); ylabel('Correct cases (%)','interpreter','tex')
else
    fprintf('\n');
    fprintf('----------------------------------------\n');
    fprintf('---------------- %s -------------------\n',index);
    fprintf('----------------------------------------\n');
    fprintf('L coverage: %i out of %i (%.2f%%)\n',sum(LCorrects),nCases,sum(LCorrects)/nCases*100);
    fprintf('NL coverage: %i out of %i (%.2f%%)\n',sum(NLCorrects),nCases,sum(NLCorrects)/nCases*100);
    fprintf('----------------------------------------\n');
    fprintf('L coverage on OR failures: %i out of %i (%.2f%%)\n',sum(LCorrects(~ORCorrects)),sum(~ORCorrects),sum(LCorrects(~ORCorrects))/sum(~ORCorrects)*100);
    fprintf('----------------------------------------\n');

    figure; hold on
    patch([0 nCases+1 nCases+1 0],[0 0 4 4],[1 1 1]);
    p(1) = bar(1:nCases,ORCorrects*4,'g');
    patch([0 nCases+1 nCases+1 0],[0 0 2 2],[1 1 1]);
    p(2) = bar(1:nCases,LCorrects*2,'r');
    patch([0 nCases+1 nCases+1 0],[0 0 1 1],[1 1 1]);
    p(3) = bar(1:nCases,NLCorrects,'y');
    yticks([]);
    title(index); xlabel('Case'); ylabel('Is correct'); axis tight; set(gcf,'position',[0,0,2000,1000]);
end
    axis tight
    legend(p,'OR','L','NL','Location','best');
    
    % Error
    output = [];
    
    errorMedian = prctile(ORError,50);
    firstQuartile = prctile(ORError,25);
    thirdQuartile = prctile(ORError,75);
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


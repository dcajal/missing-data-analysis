function [ output, p ] = computeTimeCoverage( input, errorThreshold, index )
reference = input(1,:);
or = input(2,:);
l = input(3,:);
nl = input(4,:);

nCases = length(reference);

orCorrects = zeros(nCases,numel(errorThreshold));
lCorrects = zeros(nCases,numel(errorThreshold));
nlCorrects = zeros(nCases,numel(errorThreshold));

orError = zeros(nCases,1);
lError = zeros(nCases,1);
nlError = zeros(nCases,1);
for kk=1:nCases
    orCorrects(kk,:) = abs((or{kk}.(index) - reference{kk}.(index))/reference{kk}.(index))<errorThreshold;
    lCorrects(kk,:) = abs((l{kk}.(index) - reference{kk}.(index))/reference{kk}.(index))<errorThreshold;
    nlCorrects(kk,:) = abs((nl{kk}.(index) - reference{kk}.(index))/reference{kk}.(index))<errorThreshold;
    
    % Absolute
%     orError(kk) = abs(or{kk}.(index) - reference{kk}.(index));
%     lError(kk) = abs(l{kk}.(index) - reference{kk}.(index));
%     nlError(kk) = abs(nl{kk}.(index) - reference{kk}.(index));
    
    % Relative
%     orError(kk) = abs(or{kk}.(index) - reference{kk}.(index))/reference{kk}.(index)*100;
%     lError(kk) = abs(l{kk}.(index) - reference{kk}.(index))/reference{kk}.(index)*100;
%     nlError(kk) = abs(nl{kk}.(index) - reference{kk}.(index))/reference{kk}.(index)*100;
    
    % Relative (with sign)
    orError(kk) = (or{kk}.(index) - reference{kk}.(index))/reference{kk}.(index)*100;
    lError(kk) = (l{kk}.(index) - reference{kk}.(index))/reference{kk}.(index)*100;
    nlError(kk) = (nl{kk}.(index) - reference{kk}.(index))/reference{kk}.(index)*100;
end

if numel(errorThreshold) > 1
%     figure; 
    p(1) = plot(errorThreshold*100,100*sum(orCorrects,1)/nCases,'b','Linewidth',1.5); hold on
    p(2) = plot(errorThreshold*100,100*sum(lCorrects,1)/nCases,'g','Linewidth',1.5);
    p(3) = plot(errorThreshold*100,100*sum(nlCorrects,1)/nCases,'r','Linewidth',1.5);
%     title(index); 
    xlabel('Permitted error (%)','interpreter','tex'); ylabel('Correct cases (%)','interpreter','tex')
else
    fprintf('\n');
    fprintf('----------------------------------------\n');
    fprintf('---------------- %s -------------------\n',index);
    fprintf('----------------------------------------\n');
    fprintf('OR coverage: %i out of %i (%.2f%%)\n',sum(orCorrects),nCases,sum(orCorrects)/nCases*100);
    fprintf('L coverage: %i out of %i (%.2f%%)\n',sum(lCorrects),nCases,sum(lCorrects)/nCases*100);
    fprintf('NL coverage: %i out of %i (%.2f%%)\n',sum(nlCorrects),nCases,sum(nlCorrects)/nCases*100);
    fprintf('----------------------------------------\n');
    fprintf('L coverage on OR failures: %i out of %i (%.2f%%)\n',sum(lCorrects(~orCorrects)),sum(~orCorrects),sum(lCorrects(~orCorrects))/sum(~orCorrects)*100);
    fprintf('----------------------------------------\n');

    figure; hold on
    patch([0 nCases+1 nCases+1 0],[0 0 4 4],[1 1 1]);
    p(1) = bar(1:nCases,orCorrects*4,'g');
    patch([0 nCases+1 nCases+1 0],[0 0 2 2],[1 1 1]);
    p(2) = bar(1:nCases,lCorrects*2,'r');
    patch([0 nCases+1 nCases+1 0],[0 0 1 1],[1 1 1]);
    p(3) = bar(1:nCases,nlCorrects,'y');
    yticks([]);
    title(index); xlabel('Case'); ylabel('Is correct'); axis tight; set(gcf,'position',[0,0,2000,1000]);
end
    legend(p,'OR','L','NL','Location','best');
    axis tight

    
    % Error
    output = [];
       
    errorMedian = prctile(orError,50);
    firstQuartile = prctile(orError,25);
    thirdQuartile = prctile(orError,75);
    output = [output errorMedian firstQuartile thirdQuartile];

    errorMedian = prctile(lError,50);
    firstQuartile = prctile(lError,25);
    thirdQuartile = prctile(lError,75);
    output = [output errorMedian firstQuartile thirdQuartile];

    errorMedian = prctile(nlError,50);
    firstQuartile = prctile(nlError,25);
    thirdQuartile = prctile(nlError,75);
    output = [output errorMedian firstQuartile thirdQuartile];

end


function [ p ] = computeLombCoverage( input, errorThreshold, index )
reference = input(1,:);
awNone = input(2,:);
awIncidences = input(3,:);
awIterative = input(4,:);
awIterativeNL = input(5,:);

nCases = length(reference);

noneCorrects = zeros(nCases,numel(errorThreshold));
incidencesCorrects = zeros(nCases,numel(errorThreshold));
iterativeCorrects = zeros(nCases,numel(errorThreshold));
iterativeCorrectsNL = zeros(nCases,numel(errorThreshold));
for kk=1:nCases
    noneCorrects(kk,:) = abs((awNone{kk}.(index) - reference{kk}.(index))/reference{kk}.(index))<errorThreshold;
    incidencesCorrects(kk,:) = abs((awIncidences{kk}.(index) - reference{kk}.(index))/reference{kk}.(index))<errorThreshold;
    iterativeCorrects(kk,:) = abs((awIterative{kk}.(index) - reference{kk}.(index))/reference{kk}.(index))<errorThreshold;
    iterativeCorrectsNL(kk,:) = abs((awIterativeNL{kk}.(index) - reference{kk}.(index))/reference{kk}.(index))<errorThreshold;
end

if numel(errorThreshold) > 1
    figure; hold on
    p(1) = plot(errorThreshold*100,100*sum(noneCorrects,1)/nCases,'--');
    p(2) = plot(errorThreshold*100,100*sum(incidencesCorrects,1)/nCases,'-*');
    p(3) = plot(errorThreshold*100,100*sum(iterativeCorrects,1)/nCases,'-^');
    p(4) = plot(errorThreshold*100,100*sum(iterativeCorrectsNL,1)/nCases,'-o');
    title(index); xlabel('Permited error (%)','interpreter','tex'); ylabel('Correct cases (%)','interpreter','tex')
else
    fprintf('\n');
    fprintf('----------------------------------------\n');
    fprintf('---------------- %s -------------------\n',index);
    fprintf('----------------------------------------\n');
    fprintf('No preproc coverage: %i out of %i (%.2f%%)\n',sum(noneCorrects),nCases,sum(noneCorrects)/nCases*100);
    fprintf('Incidences coverage: %i out of %i (%.2f%%)\n',sum(incidencesCorrects),nCases,sum(incidencesCorrects)/nCases*100);
    fprintf('Iterative coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrects),nCases,sum(iterativeCorrects)/nCases*100);
    fprintf('Iterative NL coverage: %i out of %i (%.2f%%)\n',sum(iterativeCorrectsNL),nCases,sum(iterativeCorrectsNL)/nCases*100);
    fprintf('----------------------------------------\n');
    fprintf('Iterative coverage on incidences failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~incidencesCorrects)),sum(~incidencesCorrects),sum(iterativeCorrects(~incidencesCorrects))/sum(~incidencesCorrects)*100);
    fprintf('Iterative coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(iterativeCorrects(~noneCorrects)),sum(~noneCorrects),sum(iterativeCorrects(~noneCorrects))/sum(~noneCorrects)*100);
    fprintf('Incidences coverage on iterative failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~iterativeCorrects)),sum(~iterativeCorrects),sum(incidencesCorrects(~iterativeCorrects))/sum(~iterativeCorrects)*100);
    fprintf('Incidences coverage on removing outliers failures: %i out of %i (%.2f%%)\n',sum(incidencesCorrects(~noneCorrects)),sum(~noneCorrects),sum(incidencesCorrects(~noneCorrects))/sum(~noneCorrects)*100);
    fprintf('----------------------------------------\n');

    figure; hold on
    patch([0 nCases+1 nCases+1 0],[0 0 4 4],[1 1 1]);
    p(1) = bar(1:nCases,noneCorrects*4,'g');
    patch([0 nCases+1 nCases+1 0],[0 0 3 3],[1 1 1]);
    p(2) = bar(1:nCases,incidencesCorrects*3,'k');
    patch([0 nCases+1 nCases+1 0],[0 0 2 2],[1 1 1]);
    p(3) = bar(1:nCases,iterativeCorrects*2,'r');
    patch([0 nCases+1 nCases+1 0],[0 0 1 1],[1 1 1]);
    p(4) = bar(1:nCases,iterativeCorrectsNL,'y');
    yticks([]);
    title(index); xlabel('Case'); ylabel('Is correct'); axis tight; set(gcf,'position',[0,0,2000,1000]);
end
    legend(p,'No preproc','Incidences','Iterative','Iterative NL','Location','bestoutside');

end


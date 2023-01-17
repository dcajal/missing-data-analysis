function fancyBoxplotAW2(input,index)
%FANCYBOXPLOTAW2 Fancy boxplots for AW dataset

for kk = 1:size(input,2)
    reference(1,kk) = input{1,kk}.(index);
    m(1,kk) = input{2,kk}.(index);
    l(1,kk) = input{3,kk}.(index);
    nl(1,kk) = input{4,kk}.(index);
end
aw = [m; l; nl]';

x1 = 1:4;

color1 = [0 166 73]/255; % Green

% figure
b1 = boxplot(reference,x1(1),'Positions',x1(1),'Colors',color1,'Widths',0.25);
hold on
b2 = boxplot(aw,x1(2:end),'Positions',x1(2:end),'Colors','b','Widths',0.25);
set(b1,'LineWidth',1.5)
set(b2,'LineWidth',1.5)

xtickangle(0)
% xtickangle(45)
set(gca,'XTick',1:4,'XTickLabel',{'Reference (ECG)', 'M', 'L', 'NL'})
% set(gca,'XTick',1:4,'XTickLabel',{'Reference (ECG)', 'OR', 'L', 'NL'})
xlim([0.5 4.5])
% ylim([0 1.1])
ylim tight

% issignificant(x1,x2,valuesRelax,valuesStress,significance);

end


function issignificant(x1,x2,valuesRelax,valuesStress,significance)
%ISSIGNIFICANT Mark as statisically significance a pair of boxplots

yLim = get(gca,'YLim');
tickHeight = 0.02*(yLim(2)-yLim(1));

hold on
for jj = 1:size(valuesRelax,2)
    if significance(jj)<0.001
        yPos = max([valuesRelax(:,jj); valuesStress(:,jj)])+1.5*tickHeight;
        line([x1(jj) x2(jj)],[yPos yPos],'Color','k','LineWidth',1.5)
        line([x1(jj) x1(jj)],[yPos yPos-tickHeight],'Color','k','LineWidth',1.5)
        line([x2(jj) x2(jj)],[yPos yPos-tickHeight],'Color','k','LineWidth',1.5)
        text(mean([x1(jj) x2(jj)])-0.09*length('**'),yPos+0.5*tickHeight,...
            '**','FontSize',20);
    elseif significance(jj)<0.05
        yPos = max([valuesRelax(:,jj); valuesStress(:,jj)])+1.5*tickHeight;
        line([x1(jj) x2(jj)],[yPos yPos],'Color','k','LineWidth',1.5)
        line([x1(jj) x1(jj)],[yPos yPos-tickHeight],'Color','k','LineWidth',1.5)
        line([x2(jj) x2(jj)],[yPos yPos-tickHeight],'Color','k','LineWidth',1.5)
        text(mean([x1(jj) x2(jj)])-0.09*length('*'),yPos+0.5*tickHeight,...
            '*','FontSize',20);
    end
end

end


function fancyBoxplotAW(reference,ro,l,nl,significance)
%FANCYBOXPLOTAW Fancy boxplots for AW dataset

valuesRelax = [reference(1:2:end)', ro(1:2:end)', l(1:2:end)', nl(1:2:end)'];
valuesStress = [reference(2:2:end)', ro(2:2:end)', l(2:2:end)', nl(2:2:end)'];

x1 = 0.825:1:4-0.175+1;
x2 = 1.175:1:4+0.175+1;
x1(2) = [];
x2(2) = [];

color1 = [0 166 73]/255; % Green

% figure
b1 = boxplot(valuesRelax,x1,'Positions',x1,'Colors',color1,'Widths',0.25);
hold on
b2 = boxplot(valuesStress,x2,'Positions',x2,'Colors','b','Widths',0.25);
set(b1,'LineWidth',1.5)
set(b2,'LineWidth',1.5)

% ylim([0 1.1])
ylim auto
xlim([0 length(x1)+2])

% xlabel('Burst duration [s]')
xtickangle(0)
% set(gca,'XTick',[1 3:size(valuesRelax,2)+1],'XTickLabel',{'Reference (ECG)', 'OR', 'L', 'NL'})
% xtickangle(45)
set(gca,'XTick',[1 3:size(valuesRelax,2)+1],'XTickLabel',{'Reference (ECG)', 'M', 'L', 'NL'})

issignificant(x1,x2,valuesRelax,valuesStress,significance);

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


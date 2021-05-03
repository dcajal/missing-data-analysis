function fancyBoxplot(values1,values2,stages,significance,mod)
%FANCYBOXPLOT Fancy boxplots

if nargin < 5
    mod = false;
end

if size(values1) ~= size(values2)
    disp('Values have different sizes')
end

x1 = 0.825:1:size(values1,2)-0.175;
x2 = 1.175:1:size(values1,2)+0.175;
if mod
    x1 = 0.825:1:size(values1,2)-0.175+1;
    x2 = 1.175:1:size(values1,2)+0.175+1;
    x1(2) = [];
    x2(2) = [];
end

color1 = [200 200 200]/255; % Grey
if mod
    color1 = [0 166 73]/255; % Green
end
    
figure('DefaultAxesFontSize',14)
b1 = boxplot(values1,x1,'Positions',x1,'Colors',color1,'Widths',0.25);
hold on
b2 = boxplot(values2,x2,'Positions',x2,'Colors','b','Widths',0.25);
set(b1,'LineWidth',1.5)
set(b2,'LineWidth',1.5)

% plot(repmat(x1',1,size(values1,1)),values1','-','Marker','.','MarkerSize',10,'Color',color1)
% plot(repmat(x2',1,size(values1,1)),values2','-','Marker','.','MarkerSize',10,'Color',color2)

% ylim([0 1.1])
ylim auto
xlim([0 length(x1)+1])
if mod
    xlim([0 length(x1)+2])
end
xlabel('Deletion probability [%]')
% xlabel('Burst duration [s]')
set(gca,'XTick',1:size(values1,2),'XTickLabel',stages)
if mod
    set(gca,'XTick',[1 3:size(values1,2)+1],'XTickLabel',stages)
end

if nargin > 3
    issignificant(x1,x2,values1,values2,significance);
end

end


function issignificant(x1,x2,values1,values2,significance)
%ISSIGNIFICANT Mark as statisically significance a pair of boxplots

yLim = get(gca,'YLim');
tickHeight = 0.02*(yLim(2)-yLim(1));

hold on
for jj = 1:size(values1,2)
    if significance(jj)<0.025
        yPos = max([values1(:,jj); values2(:,jj)])+1.5*tickHeight;
        line([x1(jj) x2(jj)],[yPos yPos],'Color','k','LineWidth',1.5)
        line([x1(jj) x1(jj)],[yPos yPos-tickHeight],'Color','k','LineWidth',1.5)
        line([x2(jj) x2(jj)],[yPos yPos-tickHeight],'Color','k','LineWidth',1.5)
        text(mean([x1(jj) x2(jj)])-0.09*length('**'),yPos+0.5*tickHeight,...
            '**','FontSize',20);
    elseif significance(jj)<0.05
        yPos = max([values1(:,jj); values2(:,jj)])+1.5*tickHeight;
        line([x1(jj) x2(jj)],[yPos yPos],'Color','k','LineWidth',1.5)
        line([x1(jj) x1(jj)],[yPos yPos-tickHeight],'Color','k','LineWidth',1.5)
        line([x2(jj) x2(jj)],[yPos yPos-tickHeight],'Color','k','LineWidth',1.5)
        text(mean([x1(jj) x2(jj)])-0.09*length('*'),yPos+0.5*tickHeight,...
            '*','FontSize',20);
    end
end

end


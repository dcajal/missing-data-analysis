function [Output] = lpp(tk,lag,removeOutliers,plotflag)
%LPP Lagged Poincare Plot

if nargin < 4
    plotflag = false;
end

% Compute RR and remove gaps via deviation from median
RR = diff(tk);
if removeOutliers
    threshold = computeThreshold(RR);
    RR(RR>threshold) = nan;
end
RR = RR*1000; % s->ms

x = RR(1:end-lag);
y = RR(1+lag:end);

% Compute LPP results
SD1 = sqrt(nanvar(x-y)/2);
SD2 = sqrt(nanvar(x+y)/2);
SD12 = SD1/SD2;
S = pi*SD1*SD2;
SDNN = sqrt((SD1^2+SD2^2)/2);

Cx = nanmean(x);
Cy = nanmean(y);
p = ([Cx Cy]-[x y]).^2;
dp = sqrt(p(:,1)+p(:,2));
Md = nanmean(dp);
Sd = nanstd(dp);

% Use a struct for output
Output.SD1 = SD1;
Output.SD2 = SD2;
Output.SD12 = SD12;
Output.S = S;
Output.SDRR = SDNN;
Output.Md = Md;
Output.Sd = Sd;

% Plot LPP
if plotflag
    figure('DefaultAxesFontSize',14);
    plot(1000*x,1000*y,'bo'); hold on;
    plot(1000*Cx,1000*Cy,'ro','LineWidth',2)
    xlabel('RR_{i} [ms]')
    ylabel('RR_{i+1} [ms]')
%     title('Lagged Poincare Plot');
    xlim([700 1200])
    ylim([700 1200])
    axis square
    line(xlim,ylim,'LineStyle','--','Color','k')
end

end


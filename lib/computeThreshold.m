function threshold = computeThreshold(rr)
    wind = 50;
    if length(rr)<wind, wind = length(rr); end
    mf = medfilt1([flipud(rr(1:wind/2)); rr; flipud(rr(end-wind/2+1:end))],wind-1);
    mf(mf>1.5) = 1.5;
    threshold = 1.5*(mf(wind/2+1:end-wind/2));
%     mf = medfilt1(rr,wind,[],'truncate');
%     threshold = 1.5*(mf);
end
function threshold = computeThreshold(rr)
    wind = 50;
    mf = medfilt1([flipud(rr(1:wind/2)); rr; flipud(rr(end-wind/2+1:end))],wind-1);
    threshold = 1.5*(mf(wind/2+1:end-wind/2));
end
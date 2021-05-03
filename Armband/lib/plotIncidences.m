function plotIncidences(signal,tk,tm,tn,fs)
    tk = 1+round(tk*fs); % samples
    % tk vector with deleted beats 
    tm = 1+round(tm*fs); % samples
    % Deleted beats 
    itm = tk(~ismember(tk,tm));
    % tk vector with deleted beats + intepolation
    tn = 1+round(tn*fs); % samples
    % Inserted beats 
    itn = tn(~ismember(tn,tm));

    ts = (0:length(signal)-1)/fs;
    figure;
    plot(ts,signal); hold on;
    plot(ts(tk+1), signal(tk+1),'ro');
    plot(ts(itn+1), signal(itn+1),'kd','markersize',10);
    plot(ts(itm+1), signal(itm+1),'kx','markersize',10);
end


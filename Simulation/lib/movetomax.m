function [newtm] = movetomax(signal,tm,fs)
% Move pulses to nearest max
% Created by Diego Cajal (2019)

    window = 0.2; % Search window [seconds]
    itm = round(tm*fs);
    
    newtm = zeros(size(itm)); % Corrected tm
    for i = 1:length(itm) % Find nearest max
        % Find nearest max at the left
        leftInterval = itm(i)+1:-1:round(itm(i)-window*fs);
        leftInterval(leftInterval<1) = [];
        [~,leftPeak] = findpeaks(signal(leftInterval),'NPeaks',1,...
            'SortStr','descend');
        if isempty(leftPeak) % No local max found
            leftPeak = window*fs+1;
        else
            leftPeak = leftPeak(1);
        end
        
        % Find nearest max at the right
        rightInterval = itm(i)-1:round(itm(i)+window*fs);
        rightInterval(rightInterval>length(signal)) = [];
        [~,rightPeak] = findpeaks(signal(rightInterval),'NPeaks',1,...
            'SortStr','descend');
        if isempty(rightPeak) % No local max found
            rightPeak = window*fs+1;
        else
            rightPeak = rightPeak(1);
        end
        
        if leftPeak < rightPeak % Choose nearest max
            newtm(i) = (leftInterval(leftPeak)-1)/fs;
        else
            newtm(i) = (rightInterval(rightPeak)-1)/fs;
        end
    end

    newtm = unique(newtm);

    figure('DefaultAxesFontSize',14)
    t = (0:length(signal)-1)/fs;
    plot(t,signal,'b'); hold on
    plot(newtm,signal(round(newtm*fs)+1),'ro')
    axis tight
end

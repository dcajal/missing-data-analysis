function [Output] = hrf(tk)
%HRF Heart Rate Fragmentation indices (from Costa 2017)

RR = diff(tk(:)');
dRR = diff(RR);
nIntervals = numel(RR);

accelerationSign = dRR;
accelerationSign(dRR>0) = 1;
accelerationSign(dRR<0) = -1;
accelerationSign(accelerationSign==0) = [];

inflectionPoints = find(diff(accelerationSign));
segmentLengths = diff(inflectionPoints);
alternationSegments = segmentLengths<2;

alternationCount = 0; % Number of intervals in alternation segments
adjacent = false;
for i = 1:(numel(alternationSegments)-1)
    if alternationSegments(i)+alternationSegments(i+1) > 1
        if adjacent % After the first 4 alternation intervals, every 1 lenght segment adds 1
            alternationCount = alternationCount+1;
        else % 4 is the minimum number of consecutive alternation intervals
            alternationCount = alternationCount+4;
            adjacent = true;
        end
    else
        adjacent = false;
    end
end

Output.pip = numel(inflectionPoints)/nIntervals*100; % Percentage of inflection points
Output.ials = 1/mean(segmentLengths); % Inverse of the average length of the acceleration/deceletration segments
Output.pss = sum(segmentLengths<3)/numel(segmentLengths)*100; % Percentage of short segments (1 or 2 beats)
Output.pas = alternationCount/nIntervals*100; % Percentage of intervals in alternation segments

end


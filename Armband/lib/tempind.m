function [Output] = tempind(tm,removeOutliers)
% TEMPIND Compute classical time domain indices
%
% Input:
%   tm: normal beat time occurrence series (in seconds)
%   removeOutliers: treat gaps as nans [Default: true]
%
% Output:
%   HRM: mean heart rate (beats/min)
%   SDNN: standard deviation of normal-to-normal (NN) intervals (ms)
%   SDSD: standard deviation of differences between adjacent NN intervals
%   RMSSD: square root of the mean squared differences of successive NN intervals
%           intervals (ms)
%   pNNx: proportion of the interval differences of successive NN intervals
%           greater that x ms with respect to all NN intervals [%]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
    removeOutliers = true;
end


% Compute RR
RR = diff(tm(:));
if removeOutliers
    threshold = computeThreshold(RR);
    RR(RR>threshold) = nan;
end
dRR = diff(RR);	% (ms)

% Compute time domain indices
HRM = nanmean(60./RR); % (beats/min)
SDNN = 1000*nanstd(RR); % (ms)
RMSSD = 1000*norm(dRR(~isnan(dRR)))/sqrt(length(dRR(~isnan(dRR)))); %(ms)
SDSD = 1000*nanstd(dRR); % (ms)
pNN50 = 100*(sum(abs(dRR)>0.050))/sum(~isnan(dRR)); % (%)    

% Use a struct for output
Output.HRM = HRM;
Output.SDNN = SDNN;
Output.SDSD = SDSD;
Output.RMSSD = RMSSD;
Output.pNN50 = pNN50;

end
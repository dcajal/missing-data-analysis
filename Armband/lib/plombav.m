function [pxx] = plombav(t,x,f,window,overlap)
% PLOMBAV    Compute Lomb periodogram averaging segments

% Inputs:
    % t: Sample times [s]
    % x: Input signal. Unevenly sampled
    % f: Frequency axis [Hz] where pxx is assessed
    % window: Segmentation window [s]
    % overlap: Overlapped seconds [s]
    
% Outputs
    % pxx: PSD estimate
    
% Created by Diego Cajal (2019)


t = t(:)';
x = x(:)';

% Compute segments
nsegments = ceil((120-window)/overlap)+1;
tsegments = zeros(2,nsegments);
tsegments(:,1) = [0 window];
for i = 2:nsegments
    tsegments(:,i) = [tsegments(1,i-1)+overlap; ...
        tsegments(1,i-1)+overlap+window];
end

% Compute Lomb
pxxAux = [];
for i = 1:nsegments
    isegment = find(t>=tsegments(1,i),1):find(t<=tsegments(2,i),1,'Last');
    if numel(isegment)>3
        pxxAux = [pxxAux plomb(x(isegment),t(isegment),f)];
    end
end

% Average segments
pxxAux = abs(pxxAux);
pxx = mean(pxxAux,2)/2;

end
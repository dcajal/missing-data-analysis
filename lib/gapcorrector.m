function tn = gapcorrector( tk,debug )
    %GAPCORRECTOR Detect and fill gaps in pulse series by interpolation
    warning('off', 'MATLAB:interp1:NaNstrip')
    
    % Threshold multipliers for upper and lower thresholds
    kupper = 0.7;
    klower = 0.6;
    
    if nargin < 2
        debug = false;
    end
    
    tk = tk(:);
    tn = tk;
    dtk = diff(tk);
    
    % Gaps are detected by deviation from the median in difference series
    threshold = computeThreshold(dtk);
    gaps = find(dtk>threshold);
    if isempty(gaps), return; end
    thresholdAtGap = threshold(gaps);
    
    % Gaps on first and last pulses are not allowed
    while gaps(1)<2
        tn(1) = [];
        dtk(1) = [];
        threshold(1) = [];
        gaps = find(dtk>threshold);
        thresholdAtGap = threshold(gaps);
    end
    while gaps(end)>numel(dtk)-1
        tn(end) = [];
        dtk(end) = [];
        threshold(end) = [];
        gaps = find(dtk>threshold);
        thresholdAtGap = threshold(gaps);
    end
     
    if debug
        f = set(gcf, 'Position', get(0, 'Screensize'));
        subplot(211);
        stem(dtk); hold on;
        stem(gaps,dtk(gaps),'r');
        hold on
        plot(threshold,'k--')
        axis tight
        ylabel('Original RR [s]')
    end
    
    nfill = 1; % Start filling with one sample
    while ~isempty(gaps)
        % In each iteration, try to fill with one more sample
        for kk = 1:length(gaps)
            auxtn = nfillgap(tn,gaps(kk),nfill);
            auxdtn = diff(auxtn);
            
            correct = auxdtn(gaps(kk):gaps(kk)+nfill)<kupper*thresholdAtGap(kk);
            limitExceeded = auxdtn(gaps(kk):gaps(kk)+nfill)<klower*thresholdAtGap(kk);
            
            if debug
                if limitExceeded
                    debugplots(auxdtn,gaps(kk),kupper*thresholdAtGap(kk),klower*thresholdAtGap(kk),nfill,false);
                else
                    debugplots(auxdtn,gaps(kk),kupper*thresholdAtGap(kk),klower*thresholdAtGap(kk),nfill,correct);
                end
            end
            
            if limitExceeded
                % Check that lower theshold is not exceeded. Use previous nfill instead
                auxtn = nfillgap(tn,gaps(kk),nfill-1);
                auxdtn = diff(auxtn);
                if debug
                    debugplots(auxdtn,gaps(kk),kupper*thresholdAtGap(kk),klower*thresholdAtGap(kk),nfill-1,true);
                end
                tn = auxtn;
                gaps = gaps+nfill-1;
            elseif correct
                % If correct number of samples, save serie
                tn = auxtn;
                gaps = gaps+nfill;
            end
        end
        
        % Compute gaps for next iteration
        dtn = diff(tn);
        threshold = computeThreshold(dtn);
        gaps = find(dtn>threshold);
        thresholdAtGap = threshold(gaps);
        nfill = nfill+1;
    end
    
    if debug
        close(f);
    end

end

function tn = nfillgap(tk,gap,nfill)
    % Interpolate gap with nfill samples    
   
    tk = [tk(1:gap); nan(nfill,1); tk(gap+1:end)];
    n = transpose(1:length(tk));
    tn = interp1(n(~isnan(tk)),tk(~isnan(tk)),n,'linear'); 
end

function debugplots(dtn,gap,upperThreshold,lowerThreshold,nfill,correct)
    subplot(212); hold off;
    stem(dtn); hold on;
    if correct
        stem(gap:gap+nfill,dtn(gap:gap+nfill),'g','LineWidth',1);
    else
        stem(gap:gap+nfill,dtn(gap:gap+nfill),'r','LineWidth',1);
    end
    axis tight
    ylabel('Corrected RR [s]')
    xlabel('Samples'); ylim([0 1.7])
    line(xlim,[upperThreshold upperThreshold],'Color','k');
    line(xlim,[lowerThreshold lowerThreshold],'Color','k');
    suptitle(sprintf('Trying with %i insertion(s)',nfill)); % sgtitle for newer matlab version
    pause;
end
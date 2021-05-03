function tn = gapcorrector( tk,debug )
    %GAPCORRECTOR Detect and fill gaps in pulse series by interpolation

    if nargin < 2
        debug = false;
    end
    
    tk = tk(:);
    tn = tk;
    dtk = diff(tk);
    
    % Gaps are detected by deviation from the median in difference series
    threshold = computeThreshold(dtk);
    gapPosition = find(dtk>threshold);
    thresholdAtGap = threshold(gapPosition);
     
    if debug
        f = figure('DefaultAxesFontSize',14);
        subplot(211);
        stem(dtk); hold on;
        stem(gapPosition,dtk(gapPosition),'r');
        hold on
        plot(threshold,'k--')
        axis tight
        ylabel('Original RR [s]')
    end
    
    nfill = 1; % Start filling with one sample
    while ~isempty(gapPosition)
        % In each iteration, try to fill with one more sample
        for kk = 1:length(gapPosition)
            auxtn = nfillgap(tn,gapPosition(kk),nfill);
            auxdtn = diff(auxtn);
            correct = auxdtn(gapPosition(kk))<0.8*thresholdAtGap(kk);
            if correct
                % If correct number of samples, save serie
                tn = auxtn;
                gapPosition = gapPosition+nfill;
            end
            
            if debug
                debugplots(auxdtn,gapPosition(kk),0.8*thresholdAtGap(kk),nfill,correct);
            end
        end
        
        % Compute gaps for next iteration
        dtn = diff(tn);
        threshold = computeThreshold(dtn);
        gapPosition = find(dtn>threshold);
        thresholdAtGap = threshold(gapPosition);
        nfill = nfill+1;
    end
    
    if debug
        close(f);
    end

end

function y = nfillgap(x,gapPosition,nfill)
    % Interpolate gapPosition with nfill samples
    
    y = [x(1:gapPosition); nan(nfill,1); x(gapPosition+1:end)];
    n = transpose(1:length(y));
    y = interp1(n(~isnan(y)),y(~isnan(y)),n,'linear');
end

function debugplots(dtn,gapPosition,med,nfill,correct)
    subplot(212); hold off;
    stem(dtn); hold on;
    if correct
        stem(gapPosition-nfill:gapPosition,dtn(gapPosition-nfill:gapPosition),'g','LineWidth',1);
    else
        stem(gapPosition:gapPosition+nfill,dtn(gapPosition:gapPosition+nfill),'r','LineWidth',1);
    end
    axis tight
    ylabel('Corrected RR [s]')
    xlabel('Samples'); ylim([0 1.8])
    line(xlim,[med med],'Color','k');
    pause;
end
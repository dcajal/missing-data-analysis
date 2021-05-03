function [z] = normalize(x)
    z = (x-nanmean(x))./nanstd(x);
end


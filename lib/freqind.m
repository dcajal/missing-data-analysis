function [ Output ] = freqind( psd, f )
HFBand = [find(f>=0.15,1) find(f>=0.4,1)];
LFBand = [find(f>=0.04,1) find(f>=0.15,1)];

Output.HF = trapz(f(HFBand(1):HFBand(2)),psd(HFBand(1):HFBand(2)));
Output.LF = trapz(f(LFBand(1):LFBand(2)),psd(LFBand(1):LFBand(2)));
Output.LFn = Output.LF/(Output.LF+Output.HF);
Output.LFHF = Output.LF/Output.HF;

end


function s_adjusted = adjustRMSbenware(s,wantedDB,BenwaredBrms)
% s_adjusted = adjustRMSbenware(s,wantedDB,BenwaredBrms)
% Scales signal s to wantedDB level.
% It supposes that in Benware, a signal with a RMS of 1 is played at
% BenwaredBrms dB (default : 94dB).

if nargin == 2,
    BenwaredBrms = 94; %dB
end

cRMS = rms(s);
wRMS = 10^((wantedDB - BenwaredBrms)/20);
s_adjusted = s .* (wRMS/cRMS);

end
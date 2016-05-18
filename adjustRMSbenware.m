function s_adjusted = adjustRMSbenware(s,wantedDB)
% s_adjusted = adjustRMSbenware(s,wantedDB,BenwaredBrms)
% Scales signal s to wantedDB level.
% It supposes that in Benware, a signal with a RMS of 1 is played at 94dB).
% 
% Reminder:
% dbSPL = 20 * log10(F/F0)
% F = F0 * 10^(dbSPL/20)

F0 = 20*(10^-6); % Pressure reference in uPascal
wRMS = F0 * 10^(wantedDB/20);
cRMS = rms(s);
s_adjusted = s .* (wRMS/cRMS);

% 
% cRMS = rms(s);
% wRMS = 10^((wantedDB - BenwaredBrms)/20);
% s_adjusted = s .* (wRMS/cRMS);

end
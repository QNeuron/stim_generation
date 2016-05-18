%% Generate pure tones for the ECoG recordings.
% 30 sec of pure tones

clear

% Pure tone parameters
% Fs           = 97656;                         % Sampling rate (Hz) (Must be at least 2 times higher than the tone frequency). 97656Hz is TDT compatible
Fs           = tdt100k;
nFreq = 20;

F            = linspace(1000,48000,nFreq);                            % Pure tone frequencies (Hz)
levels = [60 75 90]; % dB needs converion with benware --> adjustRMSbenware
benwareRefRMS = 94; % dB
nSec         = 0.1;                               % Sound duration (second)
% amplitude    = 0.25;                            % Sound amplitude (multiplicative factor)
nRep         = 20;                               % Number of repetition of the sound
interStimInt = 0.3;                             % Inter-stimulus interval (in second)
% stimName     = ['PT_' num2str(F) '_linear2'];   % Stimulus name for save file
fileDur = 30; % Wav file duration (in sec). Must not be >40sec at tdt100k for benware. ;
savePath = 'E:\\auditory-objects\\benware.stimuli\\tuningECOG_quentin\\';

fileName = 'ecogTones';

if length(dir(savePath)) > 2,
    inp = input(['The save folder (' savePath ') is not empty.\nDo you wish to continue? [N/y]'],'s');
    if isempty(inp) || strcpi(inp,'n'),
        error('Interrupted by user')
    end
end

% Create rampe
ramptype = 'linear';
dr1      = 0.025;             % Ramp duration (second)
nr1      = floor(Fs * dr1); % Number of bins covered by the ramp
dr2      = 0.025;             % Ramp duration (second)
nr2      = floor(Fs * dr2); % Number of bins covered by the ramp

switch ramptype
    case 'sinus'
        r1 = sin(linspace(0, pi/2, nr1));
        r2 = sin(linspace(0, pi/2, nr2));
    case 'linear'
        r1 = linspace(0,1,nr1);
        r2 = linspace(0,1,nr2);
    case 'gamma'
        r1 = gamma(linspace(0,2,nr1));
        r2 = linspace(0,1,nr2);
end
r = [r1, ones(1, floor((Fs*nSec) - (nr1 + nr2))), fliplr(r2)];

% Generate sound
t = linspace( 0, nSec, nSec * Fs);
c = 1;
nStim = length(levels) * length(F);
ToneMat = zeros(nStim,floor(nSec*Fs)+round(interStimInt*Fs));
ToneNFO = zeros(nStim,6); % [Frequency Level Duration Inter-stim interval]
for f = F,

% s  = amplitude * sin( 2 * pi * F * t );
s  = sin( 2 * pi * f * t );

% Apply ramp
s = s.*r;

for l = levels,
    ToneMat(c,1:floor(nSec*Fs)) = adjustRMSbenware(s,l); % Adjust dB level for benware system
    ToneNFO(c,:) = [f l nSec floor(nSec*Fs) interStimInt round(interStimInt*Fs)];
    c = c+1;
end

end

% interstimulus silence
interStimV  = zeros(nStim,round(interStimInt*Fs));
ToneMat(:,floor(nSec*Fs)+1:end) = interStimV;

% Generate 30 sec stim files
stimOrder = repmat(1:nStim,[1 nRep]);
rndind = randperm(nStim*nRep);
stimOrder = stimOrder(rndind);
totDur = (nSec + interStimInt) * nStim * nRep;
nFiles = ceil(totDur / fileDur);
nStimPerFile = fileDur/(nSec + interStimInt);

stimInfo.stimOrder = stimOrder;
stimInfo.stimInfo.info = ToneNFO;
stimInfo.stimInfo.name = {'Frequency' 'Level' 'Stim Duration (ms)' 'Stim Duration (bin)' 'Inter-stim interval (ms)' 'Inter-stim interval (bin)'};
stimInfo.fs = Fs;

indS = [1 nStimPerFile];
for i = 1:nFiles-1,
    fid = fopen(sprintf([savePath fileName '_%d.f32'],i),'w');
    for t = indS(1):indS(2),
        fwrite(fid,ToneMat(stimOrder(t),:),'float32');
    end
    localStimOrder = stimOrder(indS(1):indS(2));
    stimInfo.stimOrderF32.(sprintf([fileName '_%d'],i)) = localStimOrder;
    indS = indS + nStimPerFile;
    fclose(fid);
end

% last file
i=i+1;
fid = fopen(sprintf([savePath fileName '_%d.f32'],i),'w');
for t = indS(1):length(stimOrder),
    fwrite(fid,ToneMat(stimOrder(t),:),'float32');
end
localStimOrder = stimOrder(indS(1):length(stimOrder));
stimInfo.stimOrderF32.(sprintf([fileName '_%d'],i)) = localStimOrder;
fclose(fid);

save([savePath fileName '_stimInfo'],'stimInfo');


%%
% Add repetitions
interStimV  = zeros(1,round(interStimInt*Fs));
tmps        = [s interStimV];
s           = repmat(tmps,[1 nRep]);

% Play sound
% sound(s,Fs)

% Plot sound
figure;
t = linspace( 0, (nSec + interStimInt) * nRep, length(s));
plot(t,s);
ylim([-1.2 1.2])
set(get(gca,'xlabel'),'string','Time (Second)');
set(get(gca,'ylabel'),'string','Normalized Intensity');
title(stimName)

% Save .wav
% filename    = ['F:\matlab_CNPS\U\Users\Quentin\sounds\ABR project\' stimName] ; % Where the file will be saved
filename    = ['F:\Work\Code\Routines\Stimuli\' stimName] ; % Where the file will be saved

N           = 32;
wavwrite(s,Fs,N,filename)

%% visualising spectrogram
filename = 'F:\matlab_CNPS\U\Users\Quentin\sounds\ABR project\PT_2700_linear';
[s Fs]   = wavread(filename);

% Compute
freq                            = linspace(100,Fs/2,(Fs/2)/1000);
[~,f_SPEC,time_SPEC,SPECgram]   =spectrogram(s,1500,round(1500-Fs/1000),freq,Fs);
my                              = max(max(SPECgram));
SPECgram                        =SPECgram/my;
miny                            = 1/(10^(100/20));
SPECgram                        = max(SPECgram,miny);

%Plot
figure;
surf(time_SPEC,f_SPEC,20*log10(SPECgram),'edgecolor','none');
set(gca,'yscale','linear');axis tight; view(0,90);
colormap(1-gray.^0.7)
brighten(0.4);
grid off
set(gca,'ytick',0:2000:24000,'yticklabel',0:2:24)
set(gcf,'color',[1 1 1])



%% play around

clear

% Pure tone 1 parameter
stimName = 'test_stereo'; % Stimulus name for save file
Fs = 60000; % Sampling rate (Hz) (Must be at least 2 times higher than the tone frequency)
F = 4500; % Pure tone frequency (Hz)
nSec = 10; % Sound duration (second)
amplitude = 1; % Sound amplitude (multiplicative factor)

% Generate sound
t = linspace( 0, nSec, nSec* Fs);
s = amplitude * sin( 2 * pi * F * t );

% Pure tone 2 parameter
F2 = 9000; % Pure tone frequency (Hz)

% Generate sound
s2 = amplitude * sin( 2 * pi * F2 * t );

% sound(s,Fs)

% % Create rampe
% ramptype = 'sinus';
% dr1 = 1; % Ramp duration (second)
% nr1 = floor(Fs * dr1); % Number of bins covered by the ramp
% dr2 = 0.5; % Ramp duration (second)
% nr2 = floor(Fs * dr2); % Number of bins covered by the ramp
% 
% switch ramptype
%     case 'sinus'
%         r1 = sin(linspace(0, pi/2, nr1));
%         r2 = sin(linspace(0, pi/2, nr2));
%     case 'linear'
%         r1 = linspace(0,1,nr1);
%         r2 = linspace(0,1,nr2);
%     case 'gamma'
%         r1 = gamma(1:nr1);
%         r2 = linspace(0,1,nr2);
% end
% r = [r1, ones(1, (Fs*nSec) - (nr1 + nr2)), fliplr(r2)];

% Rampe tone 1
dr1 = 1; % Ramp duration (second)
nr1 = floor(Fs * dr1);
r1 = linspace(0,1,nr1);
dr2 = 3;
nr2 = floor(Fs * dr2);
r2 = fliplr(linspace(0,1,nr2));
dsil = 3;
nsil = floor(Fs * dsil);
sil = zeros(1,nsil);

r1 = [r1 ones(1, (Fs*nSec) - (nr1 + nr2 + nsil)) r2 sil];
r2 = fliplr(r1);

% Apply rampe
s = s.*r1;
s2 = s2.*r2;

% Play
sound([s; s2],Fs)

% Amplitude modulation

n = Fs * nSec;                        % number of samples

% set modulator
mf = 5;                            % modulator frequency (Hz)
mi = 0.5;                          % modulator index
m = (1:n) / Fs;                    % modulator data preparation
m = 1 + mi * sin(2 * pi * mf * m); % sinusoidal modulation

% amplitude modulation
s = m .* s;                        % amplitude modulation

% Play sound
sound(s,Fs)

% Plot sound
figure;
plot(t,s);
ylim(get(gca,'ylim')+0.2)
set(get(gca,'xlabel'),'string','Time (Second)');
set(get(gca,'ylabel'),'string','Normalized Intensity');

% Save .wav
filename = ['F:\matlab_CNPS\U\Users\Quentin\sounds\ABR project\' stimName] ;
N = 32;
wavwrite(s,Fs,N,filename)

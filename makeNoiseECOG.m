%% Generate Noise Burst for the ECoG recordings.
% 30 sec of Noise Burst. Every combination of rate and level will generate
% a new file, plus one file with mixed levels.

clear

% Pure tone parameters
Fs           = 97656.25;                         % Sampling rate (Hz) (Must be at least 2 times higher than the tone frequency). 97656Hz is TDT compatible
% Fs           = tdt100k;

levels = [50 60 70]; % dB needs converion with benware --> adjustRMSbenware
benwareRefRMS = 94; % dB
nSec         = 0.1;                               % Sound duration (second)
nRep         = 20;                               % Number of repetition of the sound
interStimInt = [0.3 0.1 0.5 0.75 1];                             % Inter-stimulus interval (in second)
fileDur = 30; % Wav file duration (in sec). Must not be >40sec at tdt100k for benware. ;
savePath = 'E:\\auditory-objects\\benware.stimuli\\NoiseECOG_quentin\\';

fileName = 'ecogNoise';

if length(dir(savePath)) > 2,
    inp = input(['The save folder (' savePath ') is not empty.\nDo you wish to continue? [N/y]'],'s');
    if isempty(inp)
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
% nStim = length(levels) * length(F);
% ToneMat = zeros(nStim,floor(nSec*Fs)+round(interStimInt*Fs));
% ToneNFO = zeros(nStim,6); % [Frequency Level Duration Inter-stim interval]

nFiles = length(levels) * length(interStimInt);
stimInfo.fs = Fs;
stimInfo.stimInfo.info = [];
stimInfo.stimInfo.name = {'Level' 'Stim Duration (ms)' 'Stim Duration (bin)' 'Inter-stim interval (ms)' 'Inter-stim interval (bin)'};

for i = 1:length(interStimInt)
    nStim = floor(fileDur/(nSec+interStimInt(i)));
    
    for l = levels
        
        fid = fopen(sprintf([savePath fileName '_%d.f32'],c),'w');
        for stim = 1:nStim,
            s = wgn(1,length(t),1);
            % Ramp
            s = s .* r;
            % Level adjustement
            s = adjustRMSbenware(s,l);
            % Add silence
            s = [s zeros(1,round(interStimInt(i)*Fs))];
            % Write in f32
            fwrite(fid,s,'float32');
            
        end
        fclose(fid);
        stimInfo.stimInfo.info = [stimInfo.stimInfo.info; [l nSec length(t) interStimInt(i) round(interStimInt(i)*Fs) ]];
        c = c+1;
    end
end

save([savePath fileName '_stimInfo'],'stimInfo');



%%
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


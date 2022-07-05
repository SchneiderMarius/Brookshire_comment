function [S] = Brookshire_bootstrap_01(cfg)

% Brookshire_distortion_nosampling_01 checks for the distortion that one
% gets when sampling from an AR1 distribution, due to hanning taper and
% linear detrending.
% Three signal models selcted: 'AR', 'Decay', or 'flat_LP'

if ~isfield(cfg,'ARcoef'), cfg.ARcoef = 0.5; end % AR coefficient, same as brookshire
if ~isfield(cfg, 'fs'), cfg.fs = 60 ; end % sampling rate, same as brookshire
if ~isfield(cfg, 'nSamples'), cfg.nSamples = 45; end % number of samples in trial, same as Landau
if ~isfield(cfg,'nIterations'), cfg.nIterations = 1000; end % n iterations of simulation
if ~isfield(cfg, 'scale'), cfg.scale = 0.15; end % scale factor to construct ACT, same as brookshire
if ~isfield(cfg, 'maxfreq'), cfg.maxfreq = 12; end % maximum frequency to detect
if ~isfield(cfg, 'nPermutations'), cfg.nPermutations = 5000; end
if ~isfield(cfg, 'model'), cfg.model = 'AR'; end % options: decay (in spectrum), expdecay (time domain), sinusoid, see below
if ~isfield(cfg,'nReps'), cfg.nReps = 40; end % should be even!
if ~isfield(cfg,'nSets'), cfg.nSets = 200; end % number of sets
if ~isfield(cfg,'sinmod'), cfg.sinmod = 0.2; end % number of sets
if ~isfield(cfg,'sinfreq'), cfg.sinfreq = 4; end % number of sets
if ~isfield(cfg,'freq_delta_boot'), cfg.freq_delta_boot = 1; end % tolerance for frequency deviation, 1 = +/- 1Hz
if ~isfield(cfg, 'sets'), cfg.sets = 'nonoverlapping'; end

ARcoef = cfg.ARcoef;
fs = cfg.fs; % 60 Hz.
N = cfg.nSamples; % number of data points, similar to Brookshire
nIterations = cfg.nIterations;
sc = cfg.scale; % the scale of the AR model, doesn't affect anything here - Brookshire used 0.2
maxFreq = cfg.maxfreq;
faxis = [0:N-1]/N*fs;                            %% CEM EDITS
if mod(N,2)
    faxis((N+1)/2+1:end) = faxis((N+1)/2+1:end)-fs;
else
    faxis(N/2+1:end) = faxis(N/2+1:end)-fs;
end
maxBin = find(0<=faxis&faxis<=maxFreq,1,'last'); %% CEM EDITS
nReps = cfg.nReps;

% preallocate
maxVal = NaN(1,nIterations);
indxMax = NaN(1,nIterations);
maxValBoot = NaN(nIterations, cfg.nSets); % for the halves
indxMaxBoot = NaN(nIterations,cfg.nSets); % for the halves
indxMaxTrue = NaN(1,nIterations);

% first extract statistical thresholds
pdscoreRand = NaN(1,cfg.nPermutations);
for iPerm = 1:cfg.nPermutations
    
    % equivalent to generating data according to a flat ACT
    outcome = binornd(nReps*ones(1,N),0.5*ones(1,N));
    
    % compute the ACT
    ACTest = outcome./nReps;
    
    % fft, detrend, taper
    ACTest = ACTest - nanmean(ACTest);
    ft_ACTest = fft(ACTest);
    if strcmp(cfg.stat,'max')
        pdscoreRand(iPerm) = nanmax(abs(ft_ACTest(2:maxBin)));% ./ nanmedian(abs(ftoutcome(2:end)));
    elseif strcmp(cfg.stat,'gstat')
        pdscoreRand(iPerm) = nanmax(abs(ft_ACTest(2:maxBin))) ./ nanmedian(abs(ft_ACTest(2:end)));
    end
end
srt = sort(pdscoreRand); % sort in ascending order
critval = srt(0.95*cfg.nPermutations); % 95% percentile , used for the stats

%% now generate ACTs according to an AR1
for iIter = 1:nIterations
    
    % generate a random AR1 time course
    rnd = randn(1,N);
    
    if strcmp(cfg.model, 'flat') % this just a control
        tc = ones(1,N)*0.5;
    else
        if strcmp(cfg.model, 'AR')
            tc = filter(1,[1 -ARcoef],rnd); % parameters of Brookshire, can be an option
        elseif strcmp(cfg.model,'flat_LP')
            ncomp = maxBin;
            amp = ones(1,ncomp);
            amp(1) = 1;
            X = [amp.*exp(1i*2*pi.*rand(1,ncomp)) zeros(1,N-ncomp)];
            tc = real(ifft(X));
        elseif strcmp(cfg.model, 'decay')
            ncomp = maxBin;
            amp = ones(1,ncomp)./(faxis(1:ncomp).^0.2);
            amp(1) = 1;
            X = [amp.*exp(1i*2*pi.*rand(1,ncomp)) zeros(1,N-ncomp)];
            tc = real(ifft(X));
        elseif strcmp(cfg.model, 'gauss')
            t = linspace(0,N./fs,N);
            tc = normpdf(t,t(end)/2,0.02);
        elseif strcmp(cfg.model,'chirp')
            tc = chirp(linspace(0,N./fs,N), 0, N./fs, 10);
        elseif strcmp(cfg.model, 'sin')
            t = linspace(0,N./fs,N);
            tc = sin(2*pi*cfg.sinfreq.*t+rand(1)*2*pi);
        elseif strcmp(cfg.model,'expdecay')
            t = linspace(0,N./fs,N);
            tc = exp(-t / 0.1);
        elseif strcmp(cfg.model,'decay_sin')
            ncomp = maxBin;
            amp = ones(1,ncomp)./(faxis(1:ncomp).^0.2);
            amp(1) = 1;
            amp(4) = amp(4)*(1+cfg.sinmod); % 4 Hz
            X = [amp.*exp(1i*2*pi.*rand(1,ncomp)) zeros(1,N-ncomp)];
            tc = real(ifft(X));
        elseif strcmp(cfg.model,'expdecay_sin')
            t = linspace(0,N./fs,N);
            tc = exp(-t / 0.1) + cfg.sinmod*sin(2*pi*4.*t);
        elseif strcmp(cfg.model,'gauss_sin')
            t = linspace(0,N./fs,N);
            tc = normpdf(t,t(end)/2,0.02);
            tc = tc+ cfg.sinmod*max(tc)*sin(2*pi*cfg.sinfreq.*t+rand(1)*2*pi);
        elseif strcmp(cfg.model, 'sinBrookshire')
            t = linspace(0,N./fs,N);
            tc = Brown_Noise(N);
            tc = Brookshire_rescale_timecourse(tc,sc) + ...
                sin(2*pi*6.*t+rand(1)*2*pi)*cfg.sinmod/2;
        elseif strcmp(cfg.model, 'sinAR')
            t = linspace(0,N./fs,N);
            tc = filter(1,[1 -ARcoef],rnd); % parameters of Brookshire, can be an option
            tc = Brookshire_rescale_timecourse(tc,sc) + ...
                sin(2*pi*cfg.sinfreq.*t)*cfg.sinmod;
        end
        % process the signal as in Brookshire
        if ~strcmp(cfg.model, 'sinBrookshire')
            tc = Brookshire_rescale_timecourse(tc,sc);
        end
    end
    
    % generate outcomes, compute the ACT
    outcome = zeros(nReps,N);
    for iRep = 1:nReps
        outcome(iRep,:) = binornd(ones(1,N),tc);
    end
    outcome_boot = zeros(cfg.nSets,N);
    for iBoot = 1:cfg.nSets
        for iN = 1:N % take random samples for each time point, as they are sampled in different trials
            if strcmp(cfg.sets,'nonoverlapping')
                rnd = randsample(nReps,nReps/2,false);
                x = sum(outcome(rnd,iN)) ./ (nReps/2); % mean outcome
            elseif strcmp(cfg.sets,'bootstrap')
                rnd = randsample(nReps,nReps,true);
                x = sum(outcome(rnd,iN)) ./ nReps; %(nReps/2); % mean outcome
            end
            
            outcome_boot(iBoot,iN) = x;
        end
        outcome_boot(iBoot,:) = outcome_boot(iBoot,:) - nanmean(outcome_boot(iBoot,:)); % demean already for FFT
    end
    
    % the actual outcome
    ACTest = sum(outcome)./nReps;
    
    % detrend, get rid of DC, taper, FFT
    ACTest = ACTest - nanmean(ACTest);
    ft_ACTest = fft(ACTest);
    ft_ACTboot = fft(outcome_boot,[],2);
    
    % the empirically detected maxima
    if strcmp(cfg.stat,'max')
        [maxVal(iIter),indxMax(iIter)] = nanmax(abs(ft_ACTest(2:maxBin)));
    elseif strcmp(cfg.stat,'gstat')
        [maxVal(iIter),indxMax(iIter)] = nanmax(abs(ft_ACTest(2:maxBin)) ./ nanmedian(abs(ft_ACTest(2:end))));
    end
    
    for iBoot = 1:cfg.nSets
        if strcmp(cfg.stat,'max')
            [maxValBoot(iIter,iBoot),indxMaxBoot(iIter,iBoot)] = nanmax(abs(ft_ACTboot(iBoot,2:maxBin)),[],2);
        elseif strcmp(cfg.stat,'gstat')
            [maxValBoot(iIter,iBoot),indxMaxBoot(iIter,iBoot)] = nanmax(abs(ft_ACTboot(iBoot,2:maxBin)) ./ nanmedian(abs(ft_ACTboot(iBoot,2:end))));
        end
    end
    
    % also store the actual maxima, just the spectra without sampling
    tc = tc - nanmean(tc); % demean the signal again
    ft_ACTact = fft(tc);
    
    % store
    [val,indxMaxTrue(iIter)] = nanmax(abs(ft_ACTact(2:maxBin)));
    
    % run the brookshire significance test
    [a,b,c] = armorf(ACTest-nanmean(ACTest),1,N,1);
    
    rnd = randn(1,501*N).*sqrt(b);
    tc_rnd = filter(1,[1 a],rnd); % parameters of Brookshire, can be an option
    tc_rnd = tc_rnd(N+1:end);
    tc_rnd = reshape(tc_rnd,[N 500]);
    ft = abs(fft(tc_rnd));
    mn = nanmean(abs(ft),2);
    sm = nanstd(abs(ft),[],2);
    z = (abs(ft)-mn)./sm;
    mx = nanmax(z(2:maxBin,:));
    srt = sort(mx);
    srt = srt(round(0.95*500));
    zact = (abs(ft_ACTest(:))-mn)./sm;
    zact = nanmax(zact(2:maxBin,:));
    sigBrookshire(iIter) = double(zact > srt);
    
    %  keyboard
    
end
%%
% possible criteria:
% 1) the frequency we will take is simply the empirically detected one
% - alternative is the median
% 2) how many freqs fall within some range?
% - could use all the split halves
% - or a subset of split halves
% - needs 2 parameters
[nSigFrac,fracWithin] = deal(zeros(1,cfg.nIterations));
for iIter = 1:cfg.nIterations
    
    % check how many fall within
    within = indxMaxBoot(iIter,:)>=(indxMax(iIter)-cfg.freq_delta_boot) & indxMaxBoot(iIter,:)<=(indxMax(iIter)+cfg.freq_delta_boot);
    nWithin = sum(within);
    nTot = cfg.nSets;
    fracWithin(iIter) = nWithin ./ nTot;
    
end
%

% naive detected positives
S.fracpos_naive = nansum(indxMax>1 & maxVal>critval)./ cfg.nIterations; % first frequency wouldn't count

% how many of these are really false
S.fracpos_naive_false = nansum(indxMax>1 & maxVal>critval & indxMax ~= indxMaxTrue)./ cfg.nIterations; % first frequency wouldn't count

% now with our criterion
sig_freqonly =  indxMax>1 & maxVal>critval &  fracWithin>cfg.within_freq_band_thresh;
S.fracpos_freqonly = sum(sig_freqonly) ./ cfg.nIterations;
S.fracpos_freqonly_false = sum(sig_freqonly & indxMax ~= indxMaxTrue) ./ cfg.nIterations;
S.fracpos_Brookshire = sum(sigBrookshire) ./ cfg.nIterations;
S.fracpos_Brookshire_false = sum(sigBrookshire & indxMax ~= indxMaxTrue) ./ cfg.nIterations;


function y = Brown_Noise(N)
y = cumsum(randn(1,N));

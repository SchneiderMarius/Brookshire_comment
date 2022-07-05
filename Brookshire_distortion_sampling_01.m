function [f1] = Brookshire_distortion_sampling_01(cfg)

% Brookshire_distortion_nosampling_01 checks for the distortion that one
% gets when sampling from an AR1 distribution, due to hanning taper and
% linear detrending.
% Three signal models selcted: 'AR', 'Decay', or 'flat_LP'

if nargin<1
    cfg = [];
end


if ~isfield(cfg,'ARcoef'), cfg.ARcoef = 0.5; end % AR coefficient, same as brookshire
if ~isfield(cfg, 'fs'), cfg.fs = 60 ; end % sampling rate, same as brookshire
if ~isfield(cfg, 'nsamples'), cfg.nsamples = 45; end % number of samples in trial, same as Landau
if ~isfield(cfg,'nIterations'), cfg.nIterations = 20000; end % n iterations of simulation
if ~isfield(cfg, 'scale'), cfg.scale = 0.2; end % scale factor to construct ACT, same as brookshire
if ~isfield(cfg, 'maxfreq'), cfg.maxfreq = 12; end % maximum frequency to detect
if ~isfield(cfg, 'ntrials'), cfg.ntrials = 1664; end % as in Landau
if ~isfield(cfg, 'npermutations'), cfg.npermutations = 5000; end
if ~isfield(cfg,'sinfreq'), cfg.sinfreq = 4; end % number of sets
if ~isfield(cfg, 'model'), cfg.model = 'AR'; end %


%%
lab = {'ground truth','no taper & no detrend','taper & 1st order detrend','taper & 2nd order detrend'};
col = [34, 136, 51;204, 187, 68;238, 102, 119;0,0,0]/255;
fonts   = 8;
lwidth  = 0.8;
tickl   = 0.015;

%%


ARcoef = cfg.ARcoef;
fs = cfg.fs; % 60 Hz.
N = cfg.nsamples; % number of data points, similar to Brookshire
nIterations = cfg.nIterations;
sc = cfg.scale; % the scale of the AR model, doesn't affect anything here - Brookshire used 0.2
maxFreq = cfg.maxfreq;
nTrials = cfg.ntrials;
faxis = [0:N-1]/N*fs;                            %% CEM EDITS
if mod(N,2)
    faxis((N+1)/2+1:end) = faxis((N+1)/2+1:end)-fs;
else
    faxis(N/2+1:end) = faxis(N/2+1:end)-fs;
end
maxBin = find(0<=faxis&faxis<=maxFreq,1,'last'); %% CEM EDITS
nReps = ceil(nTrials/N);


%
maxVal = NaN(3,nIterations);
indxMax = NaN(3,nIterations);
indxMaxTrue = NaN(3,nIterations);
critval = NaN(1,3);
for iProc = 1:3
    if iProc==1 % just basic
        taper = @rectwin;
        detrend_order = 0;
    elseif iProc==2 % landau
        taper = @hanning;
        detrend_order = 1;
    elseif iProc==3 % fiebelkorn
        taper = @hanning;
        detrend_order = 2;
    end
    
    % taper for computing the FFT
    w = window(taper, N);
    w = w'./sum(w);
    
    % first extract statistical thresholds
    pdscoreRand = NaN(1,cfg.npermutations);
    for iPerm = 1:cfg.npermutations
        
        % equivalent to generating data according to a flat ACT
        outcome = binornd(nReps*ones(1,N),0.5*ones(1,N));
        
        % compute the ACT
        ACTest = outcome./nReps;
        
        % fft, detrend, taper
        ACTest = detrend(ACTest,detrend_order);
        ACTest = ACTest - nanmean(ACTest);
        ft_ACTest = fft(ACTest.*w);
        pdscoreRand(iPerm) = nanmax(abs(ft_ACTest(2:maxBin)));% ./ nanmedian(abs(ftoutcome(2:end)));
    end
    srt = sort(pdscoreRand); % sort in ascending order
    critval(iProc) = srt(0.95*cfg.npermutations); % 95% percentile , used for the stats
    
    %% now generate ACTs according to an AR1
    for iIter = 1:nIterations
        
        % generate a random AR1 time course
        rnd = randn(1,N);
        
        if strcmp(cfg.model, 'flat')
            tc = ones(1,N)*0.5;
        else
            if strcmp(cfg.model, 'AR')
                tc = filter(1,[1 -ARcoef],rnd); % parameters of Brookshire, can be an option
            elseif strcmp(cfg.model,'flat_LP')
                ncomp = maxBin;
                pink=[];
                amp = ones(1,ncomp);
                amp(1) = 1;
                X = [amp.*exp(1i*2*pi.*rand(1,ncomp)) zeros(1,N-ncomp)];
                tc = real(ifft(X));
            elseif strcmp(cfg.model, 'decay')
                ncomp = maxBin;
                pink=[];
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
                tc = sin(2*pi*cfg.sinfreq.*t);
            elseif strcmp(cfg.model,'expdecay')
                t = linspace(0,N./fs,N);
                tc = exp(-t / 0.1);
            end
            
            % process the signal as in Brookshire
            tc = Brookshire_rescale_timecourse(tc,sc);
        end
        
        % generate outcomes, compute the ACT
        outcome = binornd(nReps*ones(1,N),tc);
        ACTest = outcome./nReps;
        
        % detrend, get rid of DC, taper, FFT
        ACTest = detrend(ACTest,detrend_order);
        ACTest = ACTest - nanmean(ACTest);
        ft_ACTest = fft(ACTest.*w);
        
        % the empirically detected maxima
        [maxVal(iProc,iIter),indxMax(iProc,iIter)] = nanmax(abs(ft_ACTest(2:maxBin)));
        
        % also store the actual maxima, just the spectra without sampling
        tc_detrend = detrend(tc,detrend_order); % detrend
        tc_detrend_demean = tc_detrend - nanmean(tc_detrend); % demean the signal again
        ft_ACTact = fft(tc_detrend_demean.*w);
        
        % store
        [val,indxMaxTrue(iProc,iIter)] = nanmax(abs(ft_ACTact(2:maxBin)));
    end
end

%% make histograms of the peak distribution
% plot the true peak distribution against the detected peak distribution
f1 = figure;
bins = 1:maxBin-1;

[N,B] = histc(indxMaxTrue(1,:),bins);
N = N./sum(N);
hold on
plot(faxis(bins+1),N,'Color',col(end,:),'LineWidth',lwidth)

% plot all three preprocessing ones, only the significant peaks
for iProc = 1:3
    sigIter = find(maxVal(iProc,:)>critval(iProc));
    nSig = length(sigIter);
    fracPos = nSig / nIterations;
    [N,B] = histc(indxMax(iProc,sigIter),bins);
    N = N./sum(N);
    hold on
    plot(faxis(bins+1),N,'Color',col(iProc,:),'LineWidth',lwidth)
end

% plot the ground truth
leg = legend(lab);
set(leg,'Box','off')
set(gca,'FontName', 'Arial','YTick',[0:0.5:1],'XTick',[2:2:10],'Fontsize',fonts,'Tickdir','out','TickLength',[tickl,tickl])
xlabel('Frequency [Hz]')
ylabel('Proportion')
axis square
box off
set(gcf, 'Units', 'centimeters','PaperUnits', 'centimeters', 'Renderer', 'painters')
xlim([min(faxis(bins+1)) max(faxis(bins+1))])
title(cfg.model)


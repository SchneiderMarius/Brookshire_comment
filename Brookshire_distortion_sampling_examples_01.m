function f1 = Brookshire_distortion_sampling_examples_01(cfg)

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
if ~isfield(cfg, 'model'), cfg.model = 'AR'; end %

col = [68,119,170;0, 0, 0]/255;
fonts   = 8;
lwidth  = 0.8;
tickl   = 0.015;


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
fInd = round(N/2);
maxBin = find(0<=faxis&faxis<=maxFreq,1,'last'); %% CEM EDITS
nReps = ceil(nTrials/N);

%
maxVal = NaN(3,nIterations);
indxMax = NaN(3,nIterations);
indxMaxTrue = NaN(3,nIterations);
critval = NaN(1,3);

ft_ACTestSave = [];
ACTestSave = [];
indx = [];
val = [];
for iProc = 1
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
        elseif strcmp(cfg.model,'expdecay')
            t = linspace(0,N./fs,N);
            tc = exp(-t / 0.1);
        end
        % process the signal as in Brookshire
        tc = Brookshire_rescale_timecourse(tc,sc);
        
        %% New shit
        % detrend, get rid of DC, taper, FFT
        tc2 = detrend(tc,detrend_order);
        tc2 = tc2 - nanmean(tc2);
        ft_TrueACT = fft(tc2.*w);
        %ACTfft2 = fft(tc2);
        
        
        % generate outcomes, compute the ACT
        outcome = binornd(nReps*ones(1,N),tc);
        ACTest = outcome./nReps;
        
        ACTestSave = cat(1,ACTestSave,ACTest);
        
        
        % detrend, get rid of DC, taper, FFT
        ACTest = detrend(ACTest,detrend_order);
        ACTest = ACTest - nanmean(ACTest);
        ft_ACTest = fft(ACTest.*w);
        [val(iIter),indx(iIter)] = max(abs(ft_ACTest));
        
        ft_ACTestSave = cat(1,ft_ACTestSave,abs(ft_ACTest(1:fInd)));  %% CEM EDITS
        %         ft_ACTestSave = cat(1,ft_ACTestSave,abs(ft_ACTest(1:end/2+1)));
        %ft_ACTestSave(iIter,:) = abs(ft_ACTest(1:end/2+1));
    end
end
ff = faxis;
faxis   = faxis(1:fInd);

% 4 just to select a peak around 4Hz in the example plot
id      = find(indx'==4 & max(ft_ACTestSave,[],2)>critval(iProc)); %??
[v,i]   = max(ft_ACTestSave(id(1),:));
id      = id(i);

f1 = figure;
subplot(1,2,1)
hold on
plot(linspace(0,N./fs,N), ACTestSave(id,:),'Color',col(1,:),'LineWidth',lwidth)
plot(linspace(0,N./fs,N),tc,'Color',col(2,:),'LineWidth',lwidth)
xlabel('Time [s]')
ylabel('Accuracy')
axis square
box off
set(gca,'FontName', 'Arial','YTick',[0.4:0.2:0.8],'XTick',[0:0.4:0.8],'Fontsize',fonts,'Tickdir','out','TickLength',[tickl,tickl])

subplot(1,2,2)
hold on
plot(faxis, ft_ACTestSave(id,:),'Color',col(1,:),'LineWidth',lwidth)
%plot(faxis, mean(ft_ACTestSave,1),'Color',col(2,:),'LineWidth',lwidth)
plot(faxis, abs(ft_TrueACT(1:end/2+1)),'Color',col(2,:),'LineWidth',lwidth)
plot(faxis, critval(iProc)*ones(1,length(faxis)),'k--');
xlabel('Frequency [Hz]')
ylabel('Power')
axis square
box off
leg = legend({'estimated ATC','true ATC'});
set(leg,'Box','off')
ylim([0 0.05])
set(gca,'FontName', 'Arial','YTick',[0:0.02:0.04],'XTick',[0:10:30],'Fontsize',fonts,'Tickdir','out','TickLength',[tickl,tickl])
set(gcf, 'Units', 'centimeters','PaperUnits', 'centimeters', 'Renderer', 'painters')



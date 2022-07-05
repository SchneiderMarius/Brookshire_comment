function [f] = Brookshire_distortion_nosampling_01(cfg)

% Brookshire_distortion_nosampling_01 checks for the distortion that one
% gets when sampling from an AR1 distribution, due to hanning taper and
% linear detrending.

if nargin<1
    cfg = [];
end

if ~isfield(cfg,'ARcoef'), cfg.ARcoef = 0.5; end % AR coefficient, same as brookshire
if ~isfield(cfg, 'fs'), cfg.fs = 60 ; end % sampling rate, same as brookshire
if ~isfield(cfg, 'nsamples'), cfg.nsamples = 45; end % number of samples in trial, same as brookshire
if ~isfield(cfg,'nIterations'), cfg.nIterations = 10000; end % n iterations of simulation
if ~isfield(cfg, 'scale'), cfg.scale = 0.2; end % scale factor to construct ACT, same as brookshire
if ~isfield(cfg,'sinfreq'), cfg.sinfreq = 4; end % number of sets
if ~isfield(cfg, 'maxfreq'), cfg.maxfreq = 12; end % maximum frequency to detect


%%
%col = [27,158,119;217,95,2;117,112,179]/255;
col = [0,0,0;204, 187, 68;238, 102, 119]/255;

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
faxis = [0:N-1]/N*fs;                            %% CEM EDITS
if mod(N,2)
    faxis((N+1)/2+1:end) = faxis((N+1)/2+1:end)-fs;
else
    faxis(N/2+1:end) = faxis(N/2+1:end)-fs;
end
maxBin = find(0<=faxis&faxis<=maxFreq,1,'last'); %% CEM EDITS

taperCnd = {'no taper', 'taper'};
nTaperCnd = length(taperCnd);
detrendCnd = {'no', '1st order', '2nd order'};
nDetrendCnd = length(detrendCnd);
indxMaxDistort = NaN(nTaperCnd,nDetrendCnd,nIterations);
for iIter = 1:nIterations
    
    % generate a random AR1 time course
    rnd = randn(1,N);
    tc = filter(1,[1 -ARcoef],rnd); % parameters of Brookshire, can be an option
    
    % process the signal as in Brookshire
    tc = Brookshire_rescale_timecourse(tc,sc);
    
    % the FFTs with detrending, tapering etc.
    for iTaper = 1:nTaperCnd
        
        if iTaper==1
            taper = @rectwin;
        elseif iTaper==2
            taper = @hanning;
        end
        w = window(taper, N); w = w'./sum(w);
        
        for iDetrend = 1:nDetrendCnd
            detrend_order = iDetrend - 1;
            tc_detrend = detrend(tc,detrend_order); % detrend
            tc_detrend_demean = tc_detrend - nanmean(tc_detrend); % demean the signal again
            ftest = fft(tc_detrend_demean.*w);
            
            [~,indxMaxDistort(iTaper,iDetrend,iIter)] = nanmax(abs(ftest(2:maxBin)));
        end
    end
end
%% make histograms, and then line plots
bins = 1:maxBin-1;

f = figure;

lab = {};
cmbs = [1 1; 2 2; 2 3];
for k = 1:3
    [N,B] = histc(squeeze(indxMaxDistort(cmbs(k,1),cmbs(k,2),:)),bins);
    N = N./sum(N);
    hold on
    plot(faxis(bins+1),N,'Color',col(k,:),'LineWidth',lwidth)
    str = [taperCnd{cmbs(k,1)} '_' detrendCnd{cmbs(k,2)}];
    str(strfind(str,'_')) = '-';
    if k==1
        lab{k} = 'ground truth';
    else
        lab{k} = [taperCnd{cmbs(k,1)} ' & ' detrendCnd{cmbs(k,2)} ' detrend'];
    end
end
leg = legend(lab);
set(leg,'Box','off')
set(gca,'FontName', 'Arial','Fontsize',fonts,'Tickdir','out','TickLength',[tickl,tickl])
xlabel('Frequency [Hz]')
ylabel('Proportion')
axis square
box off
set(gcf, 'Units', 'centimeters','PaperUnits', 'centimeters', 'Renderer', 'painters')
xlim([min(faxis(bins+1)) max(faxis(bins+1))])

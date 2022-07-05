function f1 = Brookshire_distortion_nosampling_examples_01(cfg)

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
if ~isfield(cfg, 'maxfreq'), cfg.maxfreq = 12; end % maximum frequency to detect


%%
col = [68,119,170;0, 0, 0]/255;

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
% faxis = linspace(0,fs,N);
% maxBin = find(faxis<=maxFreq,1,'last');
faxis = ifftshift((-N/2:N/2-1)*(fs/N));          %% CEM EDITS
maxBin = find(0<=faxis&faxis<=maxFreq,1,'last'); %% CEM EDITS

taperCnd = {'no taper'};
nTaperCnd = length(taperCnd);
detrendCnd = {'no'};
nDetrendCnd = length(detrendCnd);
indxMaxDistort = NaN(nTaperCnd,nDetrendCnd,nIterations);

ft_ACTestSave = [];
tcSave = [];

for iIter = 1:nIterations
    
    % generate a random AR1 time course
    rnd = randn(1,N);
    tc = filter(1,[1 -ARcoef],rnd); % parameters of Brookshire, can be an option
    
    % process the signal as in Brookshire
    tc = Brookshire_rescale_timecourse(tc,sc);
    
    tcSave = cat(1,tcSave,tc);
    
    % the FFTs with detrending, tapering etc.
    for iTaper = 1
        
        if iTaper==1
            taper = @rectwin;
        elseif iTaper==2
            taper = @hanning;
        end
        w = window(taper, N); w = w'./sum(w);
        
        for iDetrend = 1
            detrend_order = iDetrend - 1;
            tc_detrend = detrend(tc,detrend_order); % detrend
            tc_detrend_demean = tc_detrend - nanmean(tc_detrend); % demean the signal again
            ftest = fft(tc_detrend_demean.*w);
            
            [~,indxMaxDistort(iTaper,iDetrend,iIter)] = nanmax(abs(ftest(2:maxBin)));
            
            %         ft_ACTestSave = cat(1,ft_ACTestSave,abs(ftest(1:end/2+1)));
            ft_ACTestSave = cat(1,ft_ACTestSave,abs(ftest(1:end/2))); %% CEM EDITS
        end
    end
end


%% make histograms, and then line plots

indx = squeeze(indxMaxDistort);

% faxis   = faxis(1:end/2+1);
faxis   = faxis(1:end/2); %% CEM EDITS
id      = find(indx'==4);% & max(ft_ACTestSave,[],2)>critval(iProc));
[~,i]   = max(ft_ACTestSave(id(1),:));
id      = id(i);

f1 = figure;
subplot(1,2,1)
plot(linspace(0,N./fs,N), tcSave(id,:),'Color',col(2,:),'LineWidth',lwidth)
xlabel('Time [s]')
ylabel('Accuracy')
axis square
box off
set(gca,'FontName', 'Arial','Fontsize',fonts,'Tickdir','out','TickLength',[tickl,tickl])
set(gcf, 'Units', 'centimeters','PaperUnits', 'centimeters', 'Renderer', 'painters')

subplot(1,2,2)
hold on
plot(faxis, ft_ACTestSave(id(1),:),'Color',col(1,:),'LineWidth',lwidth)
plot(faxis, mean(ft_ACTestSave,1),'Color',col(2,:),'LineWidth',lwidth)
%plot(faxis, critval(iProc)*ones(1,length(faxis)),'k--');
xlabel('Frequency [Hz]')
ylabel('Power')
axis square
box off
leg = legend({'AR(1) realization','average'});
set(leg,'Box','off')
set(gca,'FontName', 'Arial','Fontsize',fonts,'Tickdir','out','TickLength',[tickl,tickl])

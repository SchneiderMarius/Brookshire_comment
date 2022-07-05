function [] = Brookshire_LTI_01(cfg)
% Brookshire_LTI_01 reproduces filter response and ringing

if nargin<1
    cfg = [];
end

if ~isfield(cfg, 'fs'), cfg.fs = 60 ; end % sampling rate, same as brookshire
if ~isfield(cfg, 'nsamples'), cfg.nsamples = 45; end % number of samples in trial, same as brookshire
if ~isfield(cfg, 'freqCut'), cfg.freqCut = 15; end % cut off frequency of low-pass butterworth filter
if ~isfield(cfg, 'degFilt'), cfg.degFilt = 5; end % cut off frequency of low-pass butterworth filter
if ~isfield(cfg, 'scale'), cfg.scale = 0.2; end % scale factor to construct ACT, same as brookshire
if ~isfield(cfg, 'flagSave'), cfg.flagSave = false; end % scale factor to construct ACT, same as brookshire

% fig
fonts   = 8;
figSize = 6;
tickl   = 0.015;

% params
fs = cfg.fs;
sc = cfg.scale;
N = cfg.nsamples;
flagSave = cfg.flagSave;

% filter
t = [1:N]./fs;
[b,a] = butter(cfg.degFilt,cfg.freqCut./(fs./2),'low');

% fourier
ff = [0:N-1]/N*fs;                            %% CEM EDITS
if mod(N,2)
    ff(N/2+1:end) = ff(N/2+1:end)-fs;
else
    ff((N+1)/2+1:end) = ff((N+1)/2+1:end)-fs;
end
ff = fftshift(ff);
fsel = 1 < ff & ff < 30;


% generate signals

% white
cnWhite = dsp.ColoredNoise('white','SamplesPerFrame',N*10);
rnd = cnWhite();


fA = figure,
tc = Brookshire_rescale_timecourse(rnd,sc);
sig = filter(b,a,tc);
sig = sig(N*2:N*3-1);
tc = tc(N*2:N*3-1);
plot(t,tc, 'k'), hold on
plot(t,sig, 'r'), hold on
xlabel('Time (s)')
leg = legend({'Broadband Signal','Low-Pass Filtered'});
set(leg,'Box','off')
set(gca,'FontName', 'Arial','Fontsize',fonts,'Tickdir','out','TickLength',[tickl,tickl])
axis square
ylabel('Signal Amplitude (A.U.)')
set(gcf, 'Units', 'centimeters','PaperUnits', 'centimeters', 'Renderer', 'painters', 'PaperSize', [figSize figSize], 'paperposition',[0,0,figSize,figSize])
if flagSave, saveas(fA,fullfile(cd,'FigExA.pdf'),'pdf'); end

fB = figure,
tcFourier = fftshift(fft(tc))/N;
sigFourier = fftshift(fft(sig))/N;
plot(ff(fsel), mag2db(abs(tcFourier(fsel))*2), 'k'), hold on
plot(ff(fsel), mag2db(abs(sigFourier(fsel))*2), 'r')
xlabel('Frequency (Hz)'), ylabel('Power (dB)')
xlim([min(ff(fsel)) max(ff(fsel))])
leg = legend({'Broadband Signal','Low-Pass Filtered'});
set(leg,'Box','off')
set(gca,'FontName', 'Arial','Fontsize',fonts,'Tickdir','out','TickLength',[tickl,tickl])
axis square
set(gcf, 'Units', 'centimeters','PaperUnits', 'centimeters', 'Renderer', 'painters', 'PaperSize', [figSize figSize], 'paperposition',[0,0,figSize,figSize])
if flagSave, saveas(fB,fullfile(cd,'FigExB.pdf'),'pdf'); end

fC = figure,
plot(ff(fsel), 180/pi*(angle(tcFourier(fsel))), 'k'), hold on
plot(ff(fsel), 180/pi*(angle(sigFourier(fsel))), 'r')
xlabel('Frequency (Hz)'), ylabel('Phase (Degree)')
leg = legend({'Broadband Signal','Low-Pass Filtered'});
set(leg,'Box','off')
xlim([min(ff(fsel)) max(ff(fsel))])
ylim([-180, 180])
set(gca,'FontName', 'Arial','Fontsize',fonts,'Tickdir','out','TickLength',[tickl,tickl])
axis square
set(gcf, 'Units', 'centimeters','PaperUnits', 'centimeters', 'Renderer', 'painters', 'PaperSize', [figSize figSize], 'paperposition',[0,0,figSize,figSize])
if flagSave, saveas(fC, fullfile(cd,'FigExC.pdf'),'pdf'); end

fD = figure,
tc = ones(1, N);
tc(2:end) = 0;
sig = filter(b,a,tc);
plot(t,tc, 'k'), hold on
plot(t,sig, 'r'), hold on
xlabel('Time (s)')
leg = legend({'Broadband Signal','Low-Pass Filtered'});
set(leg,'Box','off')
set(gca,'FontName', 'Arial','Fontsize',fonts,'Tickdir','out','TickLength',[tickl,tickl])
axis square
ylabel('Signal Amplitude (A.U.)')
set(gcf, 'Units', 'centimeters','PaperUnits', 'centimeters', 'Renderer', 'painters', 'PaperSize', [figSize figSize], 'paperposition',[0,0,figSize,figSize])
if flagSave, saveas(fD,fullfile(cd,'FigExD.pdf'),'pdf'); end

fE = figure,
tcFourier = fftshift(fft(tc))/N;
sigFourier = fftshift(fft(sig))/N;
plot(ff(fsel), mag2db(abs(tcFourier(fsel))*2), 'k'), hold on
plot(ff(fsel), mag2db(abs(sigFourier(fsel))*2), 'r')
xlabel('Frequency (Hz)'), ylabel('Power (dB)')
xlim([min(ff(fsel)) max(ff(fsel))])
leg = legend({'Broadband Signal','Low-Pass Filtered'});
set(leg,'Box','off')
set(gca,'FontName', 'Arial','Fontsize',fonts,'Tickdir','out','TickLength',[tickl,tickl])
axis square
set(gcf, 'Units', 'centimeters','PaperUnits', 'centimeters', 'Renderer', 'painters', 'PaperSize', [figSize figSize], 'paperposition',[0,0,figSize,figSize])
if flagSave, saveas(fE,fullfile(cd,'FigExE.pdf'),'pdf'); end


fF = figure,
plot(ff(fsel), 180/pi*(angle(tcFourier(fsel))), 'k'), hold on
plot(ff(fsel), 180/pi*(angle(sigFourier(fsel))), 'r')
xlabel('Frequency (Hz)'), ylabel('Phase (Degree)')
leg = legend({'Broadband Signal','Low-Pass Filtered'});
set(leg,'Box','off')
xlim([min(ff(fsel)) max(ff(fsel))])
ylim([-180, 180])
set(gca,'FontName', 'Arial','Fontsize',fonts,'Tickdir','out','TickLength',[tickl,tickl])
axis square
set(gcf, 'Units', 'centimeters','PaperUnits', 'centimeters', 'Renderer', 'painters', 'PaperSize', [figSize figSize], 'paperposition',[0,0,figSize,figSize])
if flagSave, saveas(fF, fullfile(cd,'FigExF.pdf'),'pdf'); end



%% Figure 1A - Example fft & Trace (no sampling)
clear  all
rng(1312);

f1 = Brookshire_distortion_nosampling_examples_01;

saveas(f1,fullfile(cd,'Fig1A.pdf'),'pdf')
close all


%% Fig 1B) - distortion & no sampling

f1 = Brookshire_distortion_nosampling_01;

set(gca,'YTick', [0:0.1:0.3])

saveas(f1,fullfile(cd,'Fig1B.pdf'),'pdf')
close all

%% Fig 1C) distortion & sampling

cfg = [];
cfg.model = 'AR';
f2 = Brookshire_distortion_sampling_01(cfg);

set(gca,'YTick', [0:0.1:0.3])

saveas(f2,fullfile(cd,'Fig1C.pdf'),'pdf')
close all
    
 
%% Fig 2A - Decay

cfg.model = 'decay';
f1 = Brookshire_distortion_sampling_examples_01(cfg);

saveas(f1,fullfile(cd,sprintf('Fig2A_%s.pdf',cfg.model)),'pdf')
close all


f1 = Brookshire_distortion_sampling_01(cfg);

saveas(f1,fullfile(cd,sprintf('Fig2B_%s.pdf',cfg.model)),'pdf')
close all

%% Fig 2B - Gauss

cfg.model = 'gauss';
f1 = Brookshire_distortion_sampling_examples_01(cfg);

saveas(f1,fullfile(cd,sprintf('Fig2C_%s.pdf',cfg.model)),'pdf')
close all


f1 = Brookshire_distortion_sampling_01(cfg);

saveas(f1,fullfile(cd,sprintf('Fig2D_%s.pdf',cfg.model)),'pdf')
close all
%% Fig 2C - Exp Decay

cfg.model = 'expdecay';
f1 = Brookshire_distortion_sampling_examples_01(cfg);

saveas(f1,fullfile(cd,sprintf('Fig2E_%s.pdf',cfg.model)),'pdf')
close all


f1 = Brookshire_distortion_sampling_01(cfg);

saveas(f1,fullfile(cd,sprintf('Fig2F_%s.pdf',cfg.model)),'pdf')
close all

%% Fig 2D

col = [102, 204, 238;170, 51, 119]/255;
fonts   = 8;
lwidth  = 0.8;
tickl   = 0.015;

models = {'AR', 'decay', 'gauss', 'expdecay' };
clear S
for iModel = 1:length(models)
  cfg.model = models{iModel};
  cfg.nIterations = 500;
  cfg.within_freq_band_thresh = 0.85; % this is the only one we really need, I use 0.9, but 0.85 could be interesting
  cfg.scale = 0.2; % suggest to show also the plot for 0.15 and 0.2 and 0.3 (for 0.3 the results are very impressive)
  cfg.nReps = 50;
  cfg.stat = 'max';
  S(iModel) = Brookshire_bootstrap_01(cfg);
end

st = []; 
for iModel = 1:length(models)
  st(iModel,2) = S(iModel).fracpos_freqonly_false; 
  st(iModel,1) = S(iModel).fracpos_naive_false;
end


f1 = figure;
b = bar(categorical(models),st);
b(1).FaceColor = col(1,:);
b(2).FaceColor = col(2,:);
box off
leg = legend({'naive rhythm detection','cross validated rhythm detection'});
set(leg,'Box','off')
set(gca,'FontName', 'Arial','YTick',[0:0.2:0.4],'Fontsize',fonts,'Tickdir','out','TickLength',[tickl,tickl])
axis square
ylabel('False positive')    
set(gcf, 'Units', 'centimeters','PaperUnits', 'centimeters', 'Renderer', 'painters')
saveas(f1,fullfile(cd,'Fig2D.pdf'),'pdf')
close all

%% Fig 2E

Modulation = [0.05 : 0.05 : 0.2];

clear Strue
cnt = 0; 
for iSc = Modulation
  cnt = cnt + 1; 
  cfg.model = 'sin'
  cfg.nIterations = 500;
  cfg.within_freq_band_thresh = 0.85; % use 0.85 or 0.9, 0.9 might be slightly too conservative
  cfg.sigfrac = 0.85;
  cfg.sig_thresh_boot = 0.85;
  cfg.scale = iSc; % suggest to show also the plot for 0.15 and 0.2 
  cfg.nReps = 50;
  cfg.stat = 'max';
  Strue(cnt) = Brookshire_bootstrap_01(cfg)
end
%
st = []; 
for ii = 1:length(Strue)
  st(ii,1) = Strue(ii).fracpos_naive; 
  st(ii,2) = Strue(ii).fracpos_freqonly;
end

f1 = figure; 
b = bar(Modulation,st);
b(1).FaceColor = col(1,:);
b(2).FaceColor = col(2,:);
box off
ylabel('Sensitivity')  
xlabel('Modulation strength')  
leg = legend({'naive rhythm detection','cross-validated rhythm detection'});
set(leg,'Box','off')
set(gca,'FontName', 'Arial','YTick',[0:0.5:1],'XTick',Modulation,'Fontsize',fonts,'Tickdir','out','TickLength',[tickl,tickl])
axis square
set(gcf, 'Units', 'centimeters','PaperUnits', 'centimeters', 'Renderer', 'painters')
saveas(f1,fullfile(cd,'Fig2E.pdf'),'pdf')
close all


%% Fig Ex 1
Brookshire_LTI_01

%% Table 1

models = {'gauss_sin','decay_sin','sin'};
clear S
scales = 0.05:0.05:0.3; 
scales = linspace(0.2, 2, 5)
for iModel = 1:3%length(models)
  iModel
  if iModel == 2
    scales = [0 1.5]
  elseif iModel ==1
    scales = [0 0.2];
  elseif iModel==3
    scales = [0 0.2];
  end
  for iSc = 2
    iSc
    cfg.model = models{iModel}; 
    cfg.nIterations = 500;
    cfg.within_freq_band_thresh = 0.85; % this is the only one we really need, I use 0.9, but 0.85 could be interesting
    cfg.scale = 0.2; % suggest to show also the plot for 0.15 and 0.2 and 0.3 (for 0.3 the results are very impressive)
    cfg.nReps = 500;
    cfg.stat = 'max';
    cfg.sinmod = scales(iSc); 
    S(iModel,iSc) = Brookshire_bootstrap_01(cfg)
  end
end

%
st = []; 
for iModel = 1:3%5%length(scales)
  st(iModel,3) = S(iModel,2).fracpos_Brookshire; 
  st(iModel,2) = S(iModel,2).fracpos_freqonly; 
  st(iModel,1) = S(iModel,2).fracpos_naive;
end

figure, 
h = bar(st)
ylim([0 1])
st = []; 
for iModel = 1:5%length(models)
  
  st(iModel,3) = S(iModel,2).fracpos_Brookshire_false; 
  st(iModel,2) = S(iModel,2).fracpos_freqonly_false; 
  st(iModel,1) = S(iModel,2).fracpos_naive_false;
end

figure, 
h = bar(st)
ylim([0 1])





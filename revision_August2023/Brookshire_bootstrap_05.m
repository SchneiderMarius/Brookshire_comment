function [S] = Brookshire_distortion_nosampling_01(cfg)

% Brookshire_distortion_nosampling_01 checks for the distortion that one
% gets when sampling from an AR1 distribution, due to hanning taper and
% linear detrending.
% Three signal models selcted: 'AR', 'Decay', or 'flat_LP'

if ~isfield(cfg,'nReps'), cfg.nReps = 40; end % should be even!
if ~isfield(cfg, 'fs'), cfg.fs = 60 ; end % sampling rate, same as brookshire
if ~isfield(cfg, 'nSamples'), cfg.nSamples = 45; end % number of samples in trial, same as Landau

% for ATC_EV
if ~isfield(cfg,'AR2par'), cfg.AR2par = [-0.9 0.2]; end
if ~isfield(cfg,'exppar'), cfg.exppar = 0.1; end % AR coefficient, same as brookshire
if ~isfield(cfg,'decaypar'), cfg.decaypar = 0.2; end % AR coefficient, same as brookshire
if ~isfield(cfg, 'model'), cfg.model = 'AR'; end % options: decay (in spectrum), expdecay (time domain), sinusoid, see below
if ~isfield(cfg, 'sinfreq'), cfg.sinfreq = 4; end
if ~isfield(cfg, 'scale'), cfg.scale = 0.2; end % scale factor to construct ACT, same as brookshire

% number of experiments
if ~isfield(cfg,'nIterations'), cfg.nIterations = 1000; end % n iterations of simulation

% for testing
if ~isfield(cfg, 'maxfreq'), cfg.maxfreq = 12; end % maximum frequency to detect
if ~isfield(cfg, 'nPermutations'), cfg.nPermutations = 5000; end
if ~isfield(cfg, 'medcv'), cfg.medcv = 'yes'; end
if ~isfield(cfg,'nMed'), cfg.nMed = 5; end

% initialize noise
rng shuffle 

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
maxBin = find(0<=faxis&faxis<=maxFreq,1,'last'); 
nReps = cfg.nReps;

% first extract statistical thresholds for shuffling in time
% this can just be done with the mean of 0.5 for all the simulations
% binomial distribution
pdscoreRand = NaN(1,cfg.nPermutations);
for iPerm = 1:cfg.nPermutations
    
    % equivalent to generating data according to a flat ACT
    outcome = binornd(nReps*ones(1,N),0.5*ones(1,N));
    
    % compute the ACT
    ACTest = outcome./nReps;
    
    % fft, detrend, taper
    ACTest = ACTest - nanmean(ACTest);
    ft_ACTest = fft(ACTest);
    pdscoreRand(iPerm) = nanmax(abs(ft_ACTest(2:maxBin)));% ./ nanmedian(abs(ftoutcome(2:end)));
end
srt = sort(pdscoreRand); % sort in ascending order
SinT.crit = srt(0.95*cfg.nPermutations); % 95% percentile , used for the stats

%
[isSigBrookshire,isSigSintPK,isSigBrookshirePK,isSigAPECS] = deal(NaN(cfg.nIterations,9));
isValid = NaN(cfg.nIterations,1);
for iIter = 1:nIterations
    iIter
    
    % generate an ATC
    str = {'decay_sin', 'expdecay_sin', 'flat_sin', 'gauss_sin', 'polynomial_sin', 'AR2_sin'}; 
    t = linspace(0,N./fs,N);      
    if  strcmp(cfg.model,str{6})
        r = zeros(1,N); r(1) = 1; 
        tc = filter(1,[1 cfg.AR2par],r); % parameters of Brookshire, can be an option     
        tc = Brookshire_rescale_timecourse(tc,sc);                                                
    elseif strcmp(cfg.model, str{3})
        tc = 0.5; 
    elseif strcmp(cfg.model, str{5})
        tc = rand.*t.^3 + rand.*t.^2 + rand.*t;  
        tc = Brookshire_rescale_timecourse(tc,sc);             
    elseif strcmp(cfg.model,str{1})
        ncomp = 10;
        amp = ones(1,ncomp)./(faxis(1:ncomp).^cfg.decaypar);
        amp(1) = 1;
                  
        X = [amp.*exp(1i*2*pi.*rand(1,ncomp)) 0.0001./(1:N-ncomp)+zeros(1,N-ncomp)];
        tc = real(ifft(X)); 
        tc = Brookshire_rescale_timecourse(tc,sc);             

    elseif strcmp(cfg.model,str{2})
        if isstr(cfg.exppar) && strcmp(cfg.exppar,'rand')
          exppar = 0.02 + rand*0.18;
        else
          exppar = cfg.exppar;
        end
        tc = exp(-t / exppar); 
        tc = Brookshire_rescale_timecourse(tc,sc);
    elseif strcmp(cfg.model,str{4})
        if isstr(cfg.gausspar) && strcmp(cfg.gausspar,'rand')
          gausspar = 0.02 + rand*0.18;
        else
          gausspar = cfg.gausspar;
        end                 
        tc = normpdf(t,t(end)/2,gausspar);
        tc = Brookshire_rescale_timecourse(tc,sc);                           
    end
     
    % set the phases equal to prevent interference ! 
    t = linspace(0,N./fs,N);           
    id = nearest(faxis,cfg.sinfreq);
    ang = angle(fft(tc)); ang = ang(id);
    ftS = angle(fft(sin(2*pi*cfg.sinfreq.*t))); angSin = ftS(id);
    tc = tc + cfg.sinmod * sin(2*pi*cfg.sinfreq.*t + (ang-angSin)); % add the sinusoid
        
    % generate outcomes, compute the ACT
    outcome = zeros(nReps,N);
    for iRep = 1:nReps
        outcome(iRep,:) = binornd(ones(1,N),tc);
    end
                 
     % accuracy time course 
    ACTest = sum(outcome)./nReps;         
    
    % get the median order
    if strcmp(cfg.medcv,'yes')
      cfg2 = []; 
      cfg2.maxBin = maxBin; 
      cfg2.meds = 5:2:9; 
      cfg2.nBoot = 50; 
      medOrd = stat_aperiodic_surrogates_medianorder(cfg2, outcome); 
    else
      medOrd = cfg.nMed; 
    end
    %%
    
    % start doing the different tests on the through to peak
    cfg2 = []; 
    cfg2.nReps = nReps; 
    cfg2.maxBin = maxBin;
    cfg2.nPerm = 500;
    cfg2.nMed = medOrd; 
    cfg2.doTaper = 0;
    cfg2.detrendord = 0;
    statAPECS = stat_test_aperiodic_surrogates(cfg2,ACTest);              
    statSINT_P2T = stat_test_SinT_P2T(cfg2,ACTest);
    statTRUE = stat_test_peakfreq(cfg2,tc);
    
    % brookshire
    cfg2 = []; 
    cfg2.nRandomizations = 500; 
    cfg2.maxBin = maxBin; 
    Brookshire = stat_test_Brookshire(cfg2,ACTest);
    
    % compute shuffling in time statistics    
    ACTest = sum(outcome)./nReps;   
    ACTest = ACTest - nanmean(ACTest);
    ft_ACTest = fft(ACTest); 
    SinT.spc = abs(ft_ACTest(2:maxBin)); 
        
    % it can occur for some gaussian parameters to find a spectral peak
    % this is not valid
    isValid(iIter) = 1; 
    if any(strcmp(cfg.model, str)) & cfg.sinmod==0 & statTRUE.indx>0
      isValid(iIter) = 0; 
    end  
    
    % if we have a sinusoid it should be there, otherwise our data is not valid.    
    if any(strcmp(cfg.model, str)) & cfg.sinmod>0 
     id = nearest(faxis(2:end), cfg.sinfreq); 
     if (statTRUE.indx-id)>0 
       isValid(iIter) = 0; 
     end
    end
      
    % collect the results
    %if isValid(iIter)==0, keyboard; end
    if any(statAPECS.indx>0)
      isSigAPECS(iIter,1) = double(any(statAPECS.pval<0.05 & statAPECS.indx==statTRUE.indx)); % TP
      isSigAPECS(iIter,2) = double(sum(statAPECS.pval>=0.05 & statAPECS.indx~=statTRUE.indx)); % CR    
      isSigAPECS(iIter,3) = double(sum(statAPECS.pval<0.05 & statAPECS.indx~=statTRUE.indx)); % FP    
      isSigAPECS(iIter,4) = double(isSigAPECS(iIter,1)==0 & statTRUE.indx>0); % false negative
      isSigAPECS(iIter,5) = sum(isSigAPECS(iIter,3))>0; % false positives trials
      isSigAPECS(iIter,6) = isSigAPECS(iIter,5)==0 & statTRUE.indx==0; % correct rejection      
      isSigAPECS(iIter,7) = double(statTRUE.indx>0); % check if there was a peak     
      if any(statAPECS.indx>0)
        isSigAPECS(iIter,8) = double(any(statAPECS.indx==statTRUE.indx)); % check if there was peak in ATC_est
      end
    else
      isSigAPECS(iIter,1:8) = 0; 
    end      
    isSigAPECS(iIter,9) = isValid(iIter);
    
    if any(statSINT_P2T.indx>0)
      issig = Brookshire.z(statSINT_P2T.indx)>Brookshire.crit; issig = issig(:)';         
      isSigBrookshirePK(iIter,1) = double(any(issig & statSINT_P2T.indx==statTRUE.indx));       
      isSigBrookshirePK(iIter,2) = double(sum(issig==0 & statSINT_P2T.indx~=statTRUE.indx)); % CR    
      isSigBrookshirePK(iIter,3) = double(sum(issig & statSINT_P2T.indx~=statTRUE.indx)); % FP    
      isSigBrookshirePK(iIter,4) = double(isSigBrookshirePK(iIter,1)==0 & statTRUE.indx>0); % FN
      isSigBrookshirePK(iIter,5) = sum(isSigBrookshirePK(iIter,3))>0;
      isSigBrookshirePK(iIter,6) = isSigBrookshirePK(iIter,5)==0 & statTRUE.indx==0; % no false positive, there was nothing       
      isSigBrookshirePK(iIter,7) = double(statTRUE.indx>0);     
      isSigBrookshirePK(iIter,8) = double(any(statSINT_P2T.indx==statTRUE.indx)); 
    else
      isSigBrookshirePK(iIter,1:8) = 0; 
    end
    isSigBrookshirePK(iIter,9) = isValid(iIter);

    if any(statSINT_P2T.indx>0)
      issig = SinT.spc(statSINT_P2T.indx)>SinT.crit; issig = issig(:)';         
      isSigSintPK(iIter,1) = double(any(issig & statSINT_P2T.indx==statTRUE.indx));       
      isSigSintPK(iIter,2) = double(sum(issig==0 & statSINT_P2T.indx~=statTRUE.indx)); % CR    
      isSigSintPK(iIter,3) = double(sum(issig & statSINT_P2T.indx~=statTRUE.indx)); % FP    
      isSigSintPK(iIter,4) = double(isSigSintPK(iIter,1)==0 & statTRUE.indx>0); % FN
      isSigSintPK(iIter,5) = sum(isSigSintPK(iIter,3))>0;
      isSigSintPK(iIter,6) = isSigSintPK(iIter,5)==0 & statTRUE.indx==0; % no false positive, there was nothing       
      isSigSintPK(iIter,7) = double(statTRUE.indx>0);     
      isSigSintPK(iIter,8) = double(any(statSINT_P2T.indx==statTRUE.indx)); 
    else
      isSigSintPK(iIter,1:6) = 0; 
    end
    isSigSintPK(iIter,9) = isValid(iIter); 
    
    issig = SinT.spc>SinT.crit; issig = issig(:)';         
    hasTruePeak = ismember(1:length(SinT.spc), statTRUE.indx);
    isSigSint(iIter,1) = double(any(issig & hasTruePeak));       
    isSigSint(iIter,2) = double(sum(issig==0 & ~hasTruePeak)); % CR    
    isSigSint(iIter,3) = double(sum(issig & ~hasTruePeak)); % FP    
    isSigSint(iIter,4) = double(isSigSint(iIter,1)==0 & statTRUE.indx>0); % FN
    isSigSint(iIter,5) = sum(isSigSint(iIter,3))>0;
    isSigSint(iIter,6) = isSigSint(iIter,5)==0 & statTRUE.indx==0; % no false positive, there was nothing       
    isSigSint(iIter,7) = double(statTRUE.indx>0);     
    isSigSint(iIter,8) = 1; 
    isSigSint(iIter,9) = isValid(iIter); 
    
    %
    issig = Brookshire.z>Brookshire.crit; issig = issig(:)';    
    hasTruePeak = ismember(1:length(Brookshire.z), statTRUE.indx);
    isSigBrookshire(iIter,1) = double(any(issig & hasTruePeak)); % TP       
    isSigBrookshire(iIter,2) = double(sum(issig==0 & ~hasTruePeak)); % CR    
    isSigBrookshire(iIter,3) = double(sum(issig & ~hasTruePeak)); % FP    
    isSigBrookshire(iIter,4) = double(isSigBrookshire(iIter,1)==0 & statTRUE.indx>0); % FN
    isSigBrookshire(iIter,5) = sum(isSigBrookshire(iIter,3))>0;
    isSigBrookshire(iIter,6) = isSigBrookshire(iIter,5)==0 & statTRUE.indx==0; % no false positive, there was nothing       
    isSigBrookshire(iIter,7) =  double(statTRUE.indx>0);
    isSigBrookshire(iIter,8) = 1; 
    isSigBrookshire(iIter,9) = isValid(iIter); 
    


end

% assemble the results
clear S
S(1,:) = sum(isSigSintPK(isValid==1,:),1); 
S(2,:) = sum(isSigBrookshirePK(isValid==1,:),1); 
S(3,:) = sum(isSigAPECS(isValid==1,:),1); 
S(4,:) = sum(isSigBrookshire(isValid==1,:),1); 
S(5,:) = sum(isSigSint(isValid==1,:),1); 

%keyboard
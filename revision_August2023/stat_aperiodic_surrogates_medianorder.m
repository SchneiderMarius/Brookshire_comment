function medOrder = stat_aperiodic_surrogates_medianorder(cfg,outcome)

% outcome is nReps x nSamples

% even number of trials
if ~isfield(cfg,'nBoot'), cfg.nBoot = 50; end
if ~isfield(cfg, 'meds'), cfg.meds = 5:2:9; end
if ~isfield(cfg, 'maxBin'), cfg.maxBin = size(outcome,2); end

nReps = size(outcome,1);

 for iBoot = 1:cfg.nBoot
   
   % take two sets of nonoverlapping trials (test and train set) 
   rnd = randsample(nReps,nReps/2,false);
   indx2 = 1:nReps;
   indx2(ismember(indx2,rnd)) = []; 
   outcome1 = outcome(rnd,:); 
   outcome2 = outcome(indx2,:);

   % get the accuracy time course
   ACTest1 = sum(outcome1)./(nReps/2);
   ACTest2 = sum(outcome2)./(nReps/2);

%       detrend, get rid of DC, taper, FFT
   ACTest1 = ACTest1 - nanmean(ACTest1);
   ft_ACTest1 = fft(ACTest1);
   ACTest2 = ACTest2 - nanmean(ACTest2);
   ft_ACTest2 = fft(ACTest2);
   cnt = 0; 
   
   % check errors relative to the median
   res = NaN(cfg.nBoot,length(cfg.meds)); 
   for iMed = cfg.meds
     cnt = cnt + 1; 
     nfft = ceil((length(ft_ACTest1)+1)/2); 
     medf = medfilt1(abs(ft_ACTest1(2:nfft)),iMed,'truncate');     
     err = abs(abs(ft_ACTest2(2:cfg.maxBin)) - medf(2:cfg.maxBin)); % remove the DC bin, 2:maxBin
     res(iBoot,cnt) = nansum(err);         
   end
   
    [val,id] = min(nanmean(res));
    medOrder = cfg.meds(id); 
 end

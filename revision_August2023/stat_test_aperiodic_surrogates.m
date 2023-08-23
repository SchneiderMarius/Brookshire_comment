function stat = stat_test_aperiodic_surrogates(cfg, ACTest)

ACTest = ACTest(:)'; 
N = length(ACTest);
mnPerformance = mean(ACTest);

maxBin = cfg.maxBin;
nPerm  = cfg.nPerm; 
nMed   = cfg.nMed;
nReps  = cfg.nReps;
dotaper = cfg.doTaper;
detrendord = cfg.detrendord;

ACTest = ACTest - nanmean(ACTest);
ft_ACTest = fft(ACTest);
nfft = ceil((length(ft_ACTest)+1)/2); 
ft = abs(ft_ACTest(2:maxBin));     

if detrendord>0
   y = detrend(ACTest,detrendord);
else
   y = ACTest;
end

if dotaper
    w = window(@hanning,N);     
    ft_ACTest = fft(y(:)'.*w(:)');
    fttest = abs(ft_ACTest(2:maxBin)); 
else
    ft_ACTest = fft(y(:)');
    fttest = abs(ft_ACTest(2:maxBin)); 
end
%keyboard
[ign,imax]=findpeaks([fttest -inf]); 
if isempty(imax)
    peak2through=0;
else
    [ign,imin]=findpeaks([-inf -fttest(1:imax(end))]);    
    peak2through = fttest(imax) - fttest(imin-1);
end
[mxpeak_est,indx_tmp] = max(peak2through);
%keyboard
if isempty(imax)
  pval = ones(1,length(fttest)); 
  stat.pval = pval;
  stat.indx = 0; 
  return;
end
 
actdata.imax= imax;
actdata.peak2through = peak2through;

%indx_max_diff = imax(indx_tmp);
mxpeak_rnd = NaN(1,nPerm);
ft_rnd = medfilt1(abs(ft_ACTest(2:nfft)),nMed,'truncate');     
for iRand = 1:nPerm
 %   keyboard
            
    spctrm = ft_rnd.*exp(1i*2*pi*rand(1,length(ft_rnd))); 
    if mod(length(ft_ACTest),2)
        tc_rnd = real(ifft([0 spctrm conj(fliplr(spctrm(1:end)))]));
    else
        tc_rnd = real(ifft([0 spctrm conj(fliplr(spctrm(1:end-1)))]));
    end            
    outcomernd = binornd(nReps*ones(1,N),tc_rnd+mnPerformance);
    ACTest_rnd = outcomernd./nReps;
    
    ACTest_rnd = ACTest_rnd - nanmean(ACTest_rnd);
    %ft_ACTest_rnd = fft(ACTest_rnd);
    %ft = abs(ft_ACTest_rnd(2:maxBin));  
    if detrendord>0
       y = detrend(ACTest_rnd,detrendord);
    else
       y = ACTest_rnd;
    end

    if dotaper
        w = window(@hanning,N);     
        ft_ACTest_rnd = fft(y(:)'.*w(:)');
        fttest = abs(ft_ACTest_rnd(2:maxBin)); 
    else
        ft_ACTest_rnd = fft(y(:)');
        fttest = abs(ft_ACTest_rnd(2:maxBin)); 
    end

  %  keyboard
    [ign,imax]=findpeaks([fttest -inf]); 
    if isempty(imax), peak2through=0;
    else
        [ign,imin]=findpeaks([-inf -fttest(1:imax(end))]);
        peak2through = fttest(imax) - fttest(imin-1);
    end

    mxpeak_rnd(iRand) = max(peak2through);
end
%keyboard
pval = zeros(1,length(actdata.imax));
indx = zeros(1,length(actdata.imax)); 
for ii = 1:length(actdata.imax)
  pval(ii) = sum((mxpeak_rnd>actdata.peak2through(ii)))/nPerm;
  indx(ii) = actdata.imax(ii);
end

stat.pval = pval;
stat.indx = indx;


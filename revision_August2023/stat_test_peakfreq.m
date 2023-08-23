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
if mxpeak_est==0
    stat.indx = 0; 
else
    indx_max_diff = imax(indx_tmp);  
    stat.indx = indx_max_diff;
end



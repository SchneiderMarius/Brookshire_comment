function Brookshire = stat_test_Brookshire(cfg,ACTest)

if ~isfield(cfg, 'nRandomizations'), cfg.nRandomizations = 500; end
if ~isfield(cfg,'maxBin'), cfg.maxBin = 15; end

ACTest = ACTest(:)'; 
N = length(ACTest);

% detrend data and compute the FFT on it
ft_ACTest_detrend = fft(detrend(ACTest,1));

% fit the AR1 model on detrended data 
[a,b,c] = armorf(detrend(ACTest,1),1,N,1);    

% generate surrogates with the same residual variance according to the AR1
% model 
rnd = randn(1,(cfg.nRandomizations+1)*N).*sqrt(b);
tc_rnd = filter(1,[1 a],rnd); % filter data with the AR1 impulse response function
tc_rnd = tc_rnd(N+1:end);
tc_rnd = reshape(tc_rnd,[N cfg.nRandomizations]);
tc_rnd = detrend(tc_rnd,1); % detrend again the data 

% get the FFT of the data and compute the zscore of the data. 
ft = abs(fft(tc_rnd));
mn = mean(abs(ft),2);
sm = std(abs(ft),[],2);
z = (abs(ft)-mn)./sm;
mx = max(z(2:cfg.maxBin,:)); % get the maximum for each randomization
srt = sort(mx); % sort the maxima
srt = srt(round(0.95*cfg.nRandomizations)); % z scores of the surrogate distribution

% actual zscore of the observed spectrum
zact = (abs(ft_ACTest_detrend(:))-mn)./sm;
z = zact(2:cfg.maxBin,:);    
[zact,ind] = max(zact(2:cfg.maxBin,:));   

Brookshire.z = z; % z-score for each frequency
Brookshire.crit = srt; % criterion for significance
     
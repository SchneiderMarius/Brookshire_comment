function tc = Brookshire_rescale_timecourse(tc, sc)

% input: 
% tc, time course
% sc, scale
% following Brookshire 2022, Nature Human Behavior

tc = tc - min(tc);
tc = (tc ./ max(tc))*sc; 
tc = tc + 0.5; 
tc = tc - nanmean(tc) + 0.5; % makes more sense


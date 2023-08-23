%% FIRST SIMULATION: DECAYING SPECTRUM, EXPLORE 4 PARAMETERS
cfgAll = {}; 
cnt = 0; 
%decaypars = [0.02 0.2]%[0.02 0.05 0.1 0.2];
decaypars = 0;
reps = [20 100 500 1000 2000];
nIter = 20;
sinmod = [0 0.02 0.04 0.06 0.08 0.1];%[0 0.01 0.02 0.04]; 
sinfreq = [4 8]; %[4 8];
nMed = 1;%
for iRep = reps
for iModel = 1:length(decaypars)
    for iIter = 1:nIter     
      for iSin = sinmod
        for iSinFreq = sinfreq
            cnt = cnt + 1; 
            cfg = []; 
            cfg.model = 'expdecay_sin';           
            cfg.exppar = 'rand';      
            cfg.nIterations = 50;
            cfg.scale = 0.2; 
            cfg.nReps = iRep;
            cfg.nPermutations = 500;
            cfg.maxfreq = 20; 
            cfg.sinmod = iSin;                   
            cfg.sinfreq = iSinFreq;
            cfg.medcv = 'yes';
        
            cfgAll{cnt} = cfg; 
          end
        end
    end
  end
end
%stat = Brookshire_bootstrap_05(cfgAll{1})
%

out = slurmfun('Brookshire_bootstrap_05',cfgAll,'partition', '8GBXS', 'waitForToolboxes', {'statistics_toolbox', 'signal_toolbox', 'image_toolbox', ...
    'curve_fitting_toolbox', 'GADS_toolbox', 'optimization_toolbox'});
save('/mnt/pns/departmentN4/MartinCode/CodeBrookshire/Brookshire_code_final/exprnd_and_sin_1','out', 'cfgAll')
%% 
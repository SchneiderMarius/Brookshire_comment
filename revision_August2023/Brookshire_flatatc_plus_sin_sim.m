%% FIRST SIMULATION: DECAYING SPECTRUM, EXPLORE 4 PARAMETERS
cfgAll = {}; 
cnt = 0; 
reps = [20 100 500 1000 2000];
nIter = 10;
sinmod = [0 0.02] %0.01 0.02 0.04]; 
sinfreq = 4 %[4 8];
nMed = 3; 
for iRep = reps
    for iModel = 1
       for iIter = 1:nIter     
          for iSin = sinmod
           for iSinFreq = sinfreq
              for iMed = 1:nMed
                cnt = cnt + 1; 
                cfg = []; 
                cfg.model = 'flat_sin';
                cfg.nIterations = 25;
                cfg.scale = 0.2; 
                cfg.nReps = iRep;
                cfg.stat = 'max';
                cfg.sinmod = iSin;           
                cfg.nPermutations = 500;
                cfg.maxfreq = 20; 
                cfg.sinfreq = iSinFreq;
                if iMed==1
                  cfg.medcv = 'yes';
                elseif iMed==2
                  cfg.medcv = 'no';
                  cfg.nMed = 5; 
                elseif iMed==2 
                  cfg.medcv = 'no';
                  cfg.nMed = 7;   
                end

                cfgAll{cnt} = cfg; 
            end
          end
        end
    end
  end
end
%
stat = Brookshire_bootstrap_04(cfgAll{1})
%%
out = slurmfun('Brookshire_bootstrap_04',cfgAll,'partition', '8GBXS', 'waitForToolboxes', {'statistics_toolbox', 'signal_toolbox', 'image_toolbox', ...
    'curve_fitting_toolbox', 'GADS_toolbox', 'optimization_toolbox'});
save('/mnt/pns/departmentN4/MartinCode/CodeBrookshire/Brookshire_code_final/flat_and_sin_1','out', 'cfgAll')
%%
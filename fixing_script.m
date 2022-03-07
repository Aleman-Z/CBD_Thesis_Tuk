addpath('C:\Users\students\Documents\Tugdual\GitHub\data_fix_0.4.2')

paths = {'F:\rat\cannabis\acutes batch 2\201\2020-04-01_12-26-15\102_CH52.continuous','F:\rat\cannabis\acutes batch 2\207\2020-04-09_11-35-16\102_CH49.continuous'  };

for ii = paths
    tic
    file = ii{1};
    fix_open_ephys_data(file);
    toc
end

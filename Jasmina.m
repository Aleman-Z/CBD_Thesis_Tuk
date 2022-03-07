%% ||||| SWRS TYPES |||||

% Script designed to run the analysis of SWRs
% Load various folders with relevants functions afterward

%Fieldtrip
addpath('C:\Users\students\Documents\Tugdual\GitHub\fieldtrip');
addpath E:\Tugdual\GitHub\fieldtrip\plotting

%CorticoHippocampal
addpath('C:\Users\students\Documents\Tugdual\GitHub\CorticoHippocampal');
addpath('C:\Users\students\Documents\Tugdual\GitHub\CorticoHippocampal\Ripple_detection');

%CBD
addpath('C:\Users\students\Documents\Tugdual\GitHub\CBD');

%downsampled file locations
addpath 'E:\rat\cannabis\PFC_all_animals'
addpath 'E:\rat\cannabis\HPC_all_animals'

%timecell
addpath E:\Tugdual\GitHub\CorticoHippocampal\General


%% Load downsampled files
%CBD animals: 
%2, 15, 16, 204, 205, 207, 209, 212, 214
%Vehicle animals:
%1, 3, 13, 14, 201, 203, 206, 208, 210, 211, 213

clear all,clc

for P=1:14
    cd('E:\rat\cannabis\PFC_all_animals')
    
    output_data_pfc = [];
    
    if P==1
        file_shal = 'PFC_201_CH7.continuous'
        file_deep = 'PFC_201_CH19.continuous'
        output_data_pfc.treat = 1; %1 = vehicle, 2 = CBD
    elseif P==2
    elseif P==3
        file_shal = 'PFC_203_CH7.continuous'
        file_deep = 'PFC_203_CH19.continuous'
        output_data_pfc.treat = 1; %1 = vehicle, 2 = CBD
    elseif P==4
        file_shal = 'PFC_204_CH7.continuous'
        file_deep = 'PFC_204_CH19.continuous'
        output_data_pfc.treat = 2; %1 = vehicle, 2 = CBD
    elseif P==5
        file_shal = 'PFC_205_CH7.continuous'
        file_deep = 'PFC_205_CH19.continuous'
        output_data_pfc.treat = 2; %1 = vehicle, 2 = CBD
    elseif P==6
        file_shal = 'PFC_206_CH7.continuous'
        file_deep = 'PFC_206_CH19.continuous'
        output_data_pfc.treat = 1; %1 = vehicle, 2 = CBD
    elseif P==7
        file_shal = 'PFC_207_CH7.continuous'
        file_deep = 'PFC_207_CH19.continuous'
        output_data_pfc.treat = 2; %1 = vehicle, 2 = CBD
    elseif P==8
        file_shal = 'PFC_208_CH7.continuous'
        file_deep = 'PFC_208_CH19.continuous'
        output_data_pfc.treat = 1; %1 = vehicle, 2 = CBD
    elseif P==9
        file_shal = 'PFC_209_CH7.continuous'
        file_deep = 'PFC_209_CH19.continuous'
        output_data_pfc.treat = 2; %1 = vehicle, 2 = CBD
    elseif P==10
        file_shal = 'PFC_210_CH7.continuous'
        file_deep = 'PFC_210_CH19.continuous'
        output_data_pfc.treat = 1; %1 = vehicle, 2 = CBD
    elseif P==11
        file_shal = 'PFC_211_CH7.continuous'
        file_deep = 'PFC_211_CH19.continuous'
        output_data_pfc.treat = 1; %1 = vehicle, 2 = CBD
    elseif P==12
        file_shal = 'PFC_212_CH7.continuous'
        file_deep = 'PFC_212_CH19.continuous'
        output_data_pfc.treat = 2; %1 = vehicle, 2 = CBD
    elseif P==13
        file_shal = 'PFC_213_CH7.continuous'
        file_deep = 'PFC_213_CH19.continuous'
        output_data_pfc.treat = 1; %1 = vehicle, 2 = CBD
    elseif P==14
        file_shal = 'PFC_214_CH7.continuous'
        file_deep = 'PFC_214_CH19.continuous'
        output_data_pfc.treat = 2; %1 = vehicle, 2 = CBD
    end
    
    %save data
    folder = 'E:\rat\cannabis\Cleaned_Data';
    save([folder 'output_data_pfc_rat' num2str(P)],'output_data_pfc')
end

for P=1:14
    cd('E:\rat\cannabis\HPC_all_animals')
    
    output_data_hpc = [];
    
    if P==1
        file_pyr = 'HPC_201_CHpyr.continuous'
        file_ab = 'HPC_201_CHab.continuous'
        file_bel = 'HPC_201_CHbel.continuous'
        output_data_hpc.treat = 1; %1 = vehicle, 2 = CBD
    elseif P==2
    elseif P==3
        file_pyr = 'HPC_203_CHpyr.continuous'
        file_ab = 'HPC_203_CHab.continuous'
        file_bel = 'HPC_203_CHbel.continuous'
        output_data_hpc.treat = 1; %1 = vehicle, 2 = CBD
    elseif P==4
        file_pyr = 'HPC_204_CHpyr.continuous'
        file_ab = 'HPC_204_CHab.continuous'
        file_bel = 'HPC_204_CHbel.continuous'
        output_data_hpc.treat = 2; %1 = vehicle, 2 = CBD
    elseif P==5
        file_pyr = 'HPC_205_CHpyr.continuous'
        file_ab = 'HPC_205_CHab.continuous'
        file_bel = 'HPC_205_CHbel.continuous'
        output_data_hpc.treat = 2; %1 = vehicle, 2 = CBD
    elseif P==6
        file_pyr = 'HPC_206_CHpyr.continuous'
        file_ab = 'HPC_206_CHab.continuous'
        file_bel = 'HPC_206_CHbel.continuous'
        output_data_hpc.treat = 1; %1 = vehicle, 2 = CBD
    elseif P==7
        file_pyr = 'HPC_207_CHpyr.continuous'
        file_ab = 'HPC_207_CHab.continuous'
        file_bel = 'HPC_207_CHbel.continuous'
        output_data_hpc.treat = 2; %1 = vehicle, 2 = CBD
    elseif P==8
        file_pyr = 'HPC_208_CHpyr.continuous'
        file_ab = 'HPC_208_CHab.continuous'
        file_bel = 'HPC_208_CHbel.continuous'
        output_data_hpc.treat = 1; %1 = vehicle, 2 = CBD
    elseif P==9
        file_pyr = 'HPC_209_CHpyr.continuous'
        file_ab = 'HPC_209_CHab.continuous'
        file_bel = 'HPC_209_CHbel.continuous'
        output_data_hpc.treat = 2; %1 = vehicle, 2 = CBD
    elseif P==10
        file_pyr = 'HPC_210_CHpyr.continuous'
        file_ab = 'HPC_210_CHab.continuous'
        file_bel = 'HPC_210_CHbel.continuous'
        output_data_hpc.treat = 1; %1 = vehicle, 2 = CBD
    elseif P==11
        file_pyr = 'HPC_211_CHpyr.continuous'
        file_ab = 'HPC_211_CHab.continuous'
        file_bel = 'HPC_211_CHbel.continuous'
        output_data_hpc.treat = 1; %1 = vehicle, 2 = CBD
    elseif P==12
        file_pyr = 'HPC_212_CHpyr.continuous'
        file_ab = 'HPC_212_CHab.continuous'
        file_bel = 'HPC_212_CHbel.continuous'
        output_data_hpc.treat = 2; %1 = vehicle, 2 = CBD
    elseif P==13
        file_pyr = 'HPC_213_CHpyr.continuous'
        file_ab = 'HPC_213_CHab.continuous'
        file_bel = 'HPC_213_CHbel.continuous'
        output_data_hpc.treat = 1; %1 = vehicle, 2 = CBD
    elseif P==14
        file_pyr = 'HPC_214_CHpyr.continuous'
        file_ab = 'HPC_214_CHab.continuous'
        file_bel = 'HPC_214_CHbel.continuous'
        output_data_hpc.treat = 2; %1 = vehicle, 2 = CBD
    end
    
    %save data
    folder = 'E:\rat\cannabis\Cleaned_Data';
    save([folder 'output_data_hpc_rat' num2str(P)],'output_data_hpc')
end

%Store data
for P=1:14
HPC = [];
PFC = [];

%load files
HPC1 = load([file_pyr '.mat']);
HPC2 = load([file_ab '.mat']); 
HPC3 = load([file_bel '.mat']);
PFC1 = load([file_shal '.mat']); 
PFC2= load([file_deep '.mat']);

%rename fields for concat to work later
mex -O RenameField.c
HPC1 = RenameField(HPC1,'HPC1','HPC');
HPC2 = RenameField(HPC2,'HPC2','HPC');
HPC3 = RenameField(HPC3,'HPC3','HPC');
PFC1 = RenameField(PFC1,'PFC1','PFC');
PFC2 = RenameField(PFC2,'PFC2','PFC');

%concatenate
HPC = [HPC1; HPC2; HPC3];
PFC = [PFC1; PFC2];

TOT = abs(HPC(2).HPC) + abs(HPC(1).HPC) + abs(HPC(3).HPC) + abs(PFC(1).PFC)+ abs(PFC(2).PFC); % Sum of amplitudes ==> To visually assess artifacts, as they will appear in every channel and add up
L = length(PFC(1).PFC');

% Remove artifacts
tr = 5670; %Visual threshold 
outliers = false(L,1);
index = 1;
while index<L
    if TOT(index)>=tr
        outliers((index-300):index+1999) = ones(2300,1); %When we detect an artifact, remove 300 datapoints prior (build up of artifact) and 2000 datapoints after (3.5 sec)
        index = index+2000;
    else
        index = index +1;
    end
end

%Filter out artifacts & replace with the mean of the channels' medians
HPC(1).HPC(outliers) = mean(median(HPC(1).HPC)); 
HPC(2).HPC(outliers) = mean(median(HPC(2).HPC)); 
HPC(3).HPC(outliers) = mean(median(HPC(3).HPC)); 
PFC(1).PFC(outliers) = mean(median(PFC(1).PFC));
PFC(2).PFC(outliers) = mean(median(PFC(2).PFC));

%save data
folder = 'E:\rat\cannabis\Cleaned_Data';
save([folder 'HPC_rat' num2str(P)],'HPC')
save([folder 'PFC_rat' num2str(P)],'PFC')
end

%% Ripples detection in pyramidal layer of HPC
for P=[1 3 4 5 6 7 8 9 10 11 12 13 14]
%load correct files   
filename=['E:\rat\cannabis\Cleaned_Data\Cleaned_DataHPC_rat' num2str(P) '.mat'];
load (filename)

%For ripple detection, we use the shallow layer of the PFC
%First we bandpass for ripple frequency range
[e,f] = butter(3, [90/300 200/300]);
ripplespyr = filtfilt(e,f,HPC(1).HPC);

yourtimevector = (1:length(ripplespyr))/600; %in seconds

% We have the start/end time of each ripple with this
thresh = mean(ripplespyr) + 5*std(ripplespyr); %threshold for ripple detection
[S, E, M] = findRipplesLisa(ripplespyr', yourtimevector, thresh, (thresh)*(1/2), [] ); %Adrian's detection

%save data
folder = 'E:\rat\cannabis\Cleaned_Data';
save([folder '\HPCpyr_rat' num2str(P)],'yourtimevector','ripplespyr','S','E','M')
end
%% reshape data to FT format
for P=[1 3 4 5 6 7 8 9 10 11 12 13 14]
%load correct files  
filename=['E:\rat\cannabis\Cleaned_Data\HPCpyr_rat' num2str(P) '.mat'];
load (filename)

%some parameter setup
ro=600; %sample freq
l=1; %not sure what this is?
freqrange=[30:1:300]; %[30:10:300]; %[0:1:30]; %frequency range of interest
toy=[-1:.01:1]; %time window 

% Adrian's function reshapes raw data to a structure that fits with FT
[p,q]=convert2fieldtrip({yourtimevector},{ripplespyr'},l,{M},{ripplespyr'},ro);
% Adrian's function that creates the time field
%store p and q for every animal and concatenate for all animals into one P
%and one Q for veh and cbd separately. similarly for spws

%then do this
timecell=create_timecell(ro,length(p));

%now put everything together in a struct
ft_datarip = [];
ft_datarip.label = { 'HPC'};
ft_datarip.fsample = 600;
ft_datarip.trial = p(1,1:end); 
ft_datarip.time = timecell(1,1:end);

%save data
folder = 'E:\rat\cannabis\Cleaned_Data\FT_data';
save([folder '\hpcdata_rip_rat' num2str(P)],'ft_datarip')
end
%% frequency analysis using FT - ripples
for P=[1 3 4 5 6 7 8 9 10 11 12 13 14]
%load correct files  
cd 'E:\rat\cannabis\Cleaned_Data\FT_data'
load(['hpcdata_rip_rat' num2str(P)],'ft_datarip')

% Compute spectrogram
cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning'; %use 'hanning' for 0-30Hz and 'dpss' for 30-300Hz
cfg.foi = [0:10:30]; %[30:10:300]; %[0:1:30]; %frequency range of interest
cfg.t_ftimwin = .02*ones(size(cfg.foi)); %use .02 for 0-30Hz TFR and .1 for 30-300Hz TFR
cfg.tapsmofrq = 10; %use 10 for 0-30Hz TFR and 30 for 30-300Hz TFR
cfg.toi = toy;
cfg.keeptrials = 'yes';
cfg.output = 'pow';

HPCripple = ft_freqanalysis(cfg, ft_datarip);

% %visualize for each rat spearately if needed
% cfg = [];
% cfg.parameter = 'powspctrm';
% cfg.colorbar = 'yes';
% cfg.title = 'PFC shallow layer ripple';
% %change color map to "fire" 
% ft_singleplotTFR(cfg, HPCripple)

%save data
folder = 'E:\rat\cannabis\Cleaned_Data\FT_data';
save([folder '\HPCripple_rat' num2str(P)],'HPCripple')
end
%% grand avrage frequency analysis for ripples
%do descriptives
for P=[1 3 4 5 6 7 8 9 10 11 12 13 14]
    load(['HPCripple_rat' num2str(P)],'HPCripple');
    P
    cfg=[]
    TFRspw = ft_freqdescriptives(cfg, HPCripple)
    grand_TFRrip{P}=TFRspw;
end

%veh - 1, 3, 6, 8, 10, 11, 13
%cbd - 4, 5, 7, 9, 12, 14
%grandaverage
cfg = [];
cfg.keepindividual = 'yes';

%grand average veh
grandavgveh = ft_freqgrandaverage(cfg,grand_TFRrip{1},grand_TFRrip{3},grand_TFRrip{6},grand_TFRrip{8},...
    grand_TFRrip{10},grand_TFRrip{11},grand_TFRrip{13});

%grand average cbd
grandavgcbd = ft_freqgrandaverage(cfg,grand_TFRrip{4},grand_TFRrip{5},grand_TFRrip{7},grand_TFRrip{9},...
    grand_TFRrip{12},grand_TFRrip{14});

clear grand_TFRrip %free up some memory

%visualize vehicle grand average
cfg = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.title = 'VEH HPC pyramidal layer ripples 0-30Hz';
ft_colormap(hot)
ft_singleplotTFR(cfg, grandavgveh)

%visualize cbd grand average
cfg = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.title = 'CBD HPC pyramidal layer ripples 0-30Hz';
ft_colormap(hot)
ft_singleplotTFR(cfg, grandavgcbd)

%% group analysis without grand average
p_combined_veh =[];
p_combined_cbd =[];

%create concatenated trials variable
for P=[1, 3, 6, 8, 10, 11, 13] %vehicle rats
cd 'E:\rat\cannabis\Cleaned_Data\FT_data'
load(['hpcdata_rip_rat' num2str(P)],'ft_datarip')

%take the .time from each file and put into one big variable for veh and cbd rats
p_combined_veh = [p_combined_veh ft_datarip.trial];
end

%create combined timecell
timecell_combined_veh=create_timecell(ro,length(p_combined_veh));

%new concatenated struct for VEH rats
ft_vehrip = [];
ft_vehrip.label = { 'HPC_VEH'};
ft_vehrip.fsample = 600;
ft_vehrip.trial = p_combined_veh; 
ft_vehrip.time = timecell_combined_veh;

%Compute spectrogram for VEH rats
cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning'; %use 'Hanning for 0-30Hz and ''dpss' for 30-300Hz
cfg.foi = [0:10:30]; %[30:10:300]; %[0:1:30]; %frequency range of interest
cfg.t_ftimwin = .02*ones(size(cfg.foi)); %use .02 for 0-30Hz TFR and .1 for 30-300Hz TFR
cfg.tapsmofrq = 10; %use 10 for 0-30Hz TFR and 30 for 30-300Hz TFR
cfg.toi = toy;
cfg.keeptrials = 'yes';
cfg.output = 'pow';

HPCripVEH = ft_freqanalysis(cfg, ft_vehrip);

%visualize 
cfg = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.title = 'VEH HPC pyramidal layer ripples 0-30Hz';
ft_colormap(hot)
ft_singleplotTFR(cfg, HPCripVEH)


%create concatenated trials variable
for P=[4, 5, 7, 9, 12, 14] %CBD rats
cd 'E:\rat\cannabis\Cleaned_Data\FT_data'
load(['hpcdata_rip_rat' num2str(P)],'ft_datarip')

%take the .time from each file and put into one big variable for veh and cbd rats
p_combined_cbd = [p_combined_cbd ft_datarip.trial];
end

%create combined timecell
timecell_combined_cbd=create_timecell(ro,length(p_combined_cbd));

%new concatenated struct for CBD rats
ft_cbdrip = [];
ft_cbdrip.label = { 'HPC_CBD'};
ft_cbdrip.fsample = 600;
ft_cbdrip.trial = p_combined_cbd; 
ft_cbdrip.time = timecell_combined_cbd;

% Compute spectrogram for CBD rats
cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'dpss'; %use 'Hanning for 0-30Hz and ''dpss' for 30-300Hz
cfg.foi = [30:10:300]; %[30:10:300]; %[0:1:30]; %frequency range of interest
cfg.t_ftimwin = .1*ones(size(cfg.foi)); %use .02 for 0-30Hz TFR and .1 for 30-300Hz TFR
cfg.tapsmofrq = 30; %use 10 for 0-30Hz TFR and 30 for 30-300Hz TFR
cfg.toi = toy;
cfg.keeptrials = 'yes';
cfg.output = 'pow';

HPCripCBD = ft_freqanalysis(cfg, ft_cbdrip);

%visualize 
cfg = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.title = 'CBD HPC pyramidal layer ripples 30-300Hz';
ft_colormap(hot)
ft_singleplotTFR(cfg, HPCripCBD)

%% Sharp wave detection in HPC3
for P=[1 3 4 5 6 7 8 9 10 11 12 13 14]
%load correct files   
filename=['E:\rat\cannabis\Cleaned_Data\Cleaned_DataHPC_rat' num2str(P) '.mat'];
load (filename)

% For the detection of SW (sharp waves) we use the deep? layer of PFC
%Band pass filter to better see the sharp waves
[c,d] = butter(3, [1/300 20/300]);
filtbpl = filtfilt(c,d,HPC(2).HPC);
spw = double(filtbpl <= mean(filtbpl)-5*std(filtbpl)); % Detect sharp waves as signal lower than -5*SD
dspw = abs(diff(spw)); % Use absolute derivative for practical purposes: every sharp wave will hence start with a 1 and end with a 1

%Get the local peak of each sharp wave 
spwpeak = nan(L,1);
index = 1;
Indmin = [];
while index < L
    if dspw(index) == 1
        index2 = index;
        index = index+1;
        while dspw(index) == 0
            index = index + 1;
        end
        locmin = min(HPC(2).HPC(index2+1):(index+1));
        indmin = find(HPC(2).HPC == locmin);
        Indmin = [Indmin indmin];
        spwpeak(indmin) = locmin;
        index = index+1;
    else
        index = index+1;
    end
end

%save data
folder = 'E:\rat\cannabis\Cleaned_Data';
save([folder '\HPCdeep_rat' num2str(P)],'filtbpl','spwpeak','Indmin')
end
%% Reshape data to FT format
for P=[1 3 4 5 6 7 8 9 10 11 12 13 14]
%load correct files  
filename=['E:\rat\cannabis\Cleaned_Data\HPCdeep_rat' num2str(P) '.mat'];
load (filename)

%some parameter setup
ro=600; %sample freq
l=1; %not sure what this is?
freqrange=[30:1:300]; %[30:10:300]; %frequency range of interest
toy=[-1:.01:1]; %time window
timevectorsw = (1:length(spwpeak))/600; %in seconds
spwpk = Indmin./600; %M values, i.e. peaks of sw
swpfc = filtbpl./600;

% Adrian's function reshapes raw data to a structure that fits with FT
[p,q]=convert2fieldtrip({yourtimevector},{swpfc'},l,{spwpk'},{swpfc'},ro); 
% Adrian's function that creates the time field
timecell=create_timecell(ro,length(p));

%now put everything together in a struct
ft_dataspw = [];
ft_dataspw.label = { 'HPC'};
ft_dataspw.fsample = 600;
ft_dataspw.trial = p(1,1:end); 
ft_dataspw.time = timecell(1,1:end);

%save data
folder = 'E:\rat\cannabis\Cleaned_Data\FT_data';
save([folder '\hpcdata_spw_rat' num2str(P)],'ft_dataspw')
end
%% frequency analysis using FT - sharp waves
for P=[1 3 4 5 6 7 8 9 10 11 12 13 14]
%load correct files  
cd 'E:\rat\cannabis\Cleaned_Data\FT_data'
load(['hpcdata_spw_rat' num2str(P)],'ft_dataspw')

% Compute spectrogram
cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'dpss'; %use 'Hanning for 0-30Hz and ''dpss' for 30-300Hz
cfg.foi = [30:10:300]; %[30:10:300]; %[0:1:30]; %frequency range of interest
cfg.t_ftimwin = .1*ones(size(cfg.foi)); %use .02 for 0-30Hz TFR and .1 for 30-300Hz TFR
cfg.tapsmofrq = 30; %use 10 for 0-30Hz TFR and 30 for 30-300Hz TFR
cfg.toi = toy;
cfg.keeptrials = 'yes';
cfg.output = 'pow';

HPCspw = ft_freqanalysis(cfg, ft_dataspw);

% %visualize separately for each rat if needed
% cfg = [];
% cfg.parameter = 'powspctrm';
% cfg.colorbar = 'yes';
% cfg.title = 'PFC sharp waves';
% %change heat map
% ft_singleplotTFR(cfg, PFCspw)

%save data
folder = 'E:\rat\cannabis\Cleaned_Data\FT_data';
save([folder '\HPCspw_rat' num2str(P)],'HPCspw')
end
%% grand average frequency analysis for spw
%do descriptives
for P=[1 3 4 5 6 7 8 9 10 11 12 13 14]
    load(['HPCspw_rat' num2str(P)],'HPCrspw');
    P
    cfg=[]
    TFRspw = ft_freqdescriptives(cfg, HPCspw)
    grand_TFRspw{P}=TFRspw;
end

%veh - 1, 3, 6, 8, 10, 11, 13
%cbd - 4, 5, 7, 9, 12, 14
%grandaverage
cfg = [];
cfg.keepindividual = 'yes';

%grand average veh
grandavgveh = ft_freqgrandaverage(cfg,grand_TFRspw{1},grand_TFRspw{3},grand_TFRspw{6},grand_TFRspw{8},...
    grand_TFRspw{10},grand_TFRspw{11},grand_TFRspw{13});

%grand average cbd
grandavgcbd = ft_freqgrandaverage(cfg,grand_TFRspw{4},grand_TFRspw{5},grand_TFRspw{7},grand_TFRspw{9},...
    grand_TFRspw{12},grand_TFRspw{14});

clear grand_TFRspw %free up some memory

%visualize vehicle grand average
cfg = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.title = 'VEH HPC below pyramidal layer sharp wave 30-300Hz';
ft_colormap(hot)
ft_singleplotTFR(cfg, grandavgveh)

%visualize cbd grand average
cfg = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.title = 'CBD HPC below pyramidal layer sharp wave 30-300Hz';
ft_colormap(hot)
ft_singleplotTFR(cfg, grandavgcbd)

%% group analysis without grand average
p_combined_veh =[];
p_combined_cbd =[];

%create concatenated trials variable
for P=[1, 3, 6, 8, 10, 11, 13] %vehicle rats
cd 'E:\rat\cannabis\Cleaned_Data\FT_data'
load(['hpcdata_spw_rat' num2str(P)],'ft_dataspw')

%take the .time from each file and put into one big variable for veh and cbd rats
p_combined_veh = [p_combined_veh ft_dataspw.trial];
end

%create combined timecell
timecell_combined_veh=create_timecell(ro,length(p_combined_veh));

%new concatenated struct for VEH rats
ft_vehspw = [];
ft_vehspw.label = { 'HPC_VEH'};
ft_vehspw.fsample = 600;
ft_vehspw.trial = p_combined_veh; 
ft_vehspw.time = timecell_combined_veh;

%Compute spectrogram for VEH rats
cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'dpss'; %use 'Hanning for 0-30Hz and ''dpss' for 30-300Hz
cfg.foi = [30:10:300]; %[30:10:300]; %[0:1:30]; %frequency range of interest
cfg.t_ftimwin = .1*ones(size(cfg.foi)); %use .02 for 0-30Hz TFR and .1 for 30-300Hz TFR
cfg.tapsmofrq = 30; %use 10 for 0-30Hz TFR and 30 for 30-300Hz TFR
cfg.toi = toy;
cfg.keeptrials = 'yes';
cfg.output = 'pow';

HPCspwVEH = ft_freqanalysis(cfg, ft_vehspw);

%visualize 
cfg = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.title = 'VEH HPC below pyramidal layer sharp wave 30-300Hz';
ft_colormap(hot)
ft_singleplotTFR(cfg, HPCspwVEH)


%create concatenated trials variable
for P=[4, 5, 7, 9, 12, 14] %CBD rats
cd 'E:\rat\cannabis\Cleaned_Data\FT_data'
load(['hpcdata_spw_rat' num2str(P)],'ft_dataspw')

%take the .time from each file and put into one big variable for veh and cbd rats
p_combined_cbd = [p_combined_cbd ft_dataspw.trial];
end

%create combined timecell
timecell_combined_cbd=create_timecell(ro,length(p_combined_cbd));

%new concatenated struct for CBD rats
ft_cbdspw = [];
ft_cbdspw.label = {'HPC_CBD'};
ft_cbdspw.fsample = 600;
ft_cbdspw.trial = p_combined_cbd; 
ft_cbdspw.time = timecell_combined_cbd;

%Compute spectrogram for CBD rats
cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'dpss'; %use 'Hanning for 0-30Hz and ''dpss' for 30-300Hz
cfg.foi = [30:10:300]; %[30:10:300]; %[0:1:30]; %frequency range of interest
cfg.t_ftimwin = .1*ones(size(cfg.foi)); %use .02 for 0-30Hz TFR and .1 for 30-300Hz TFR
cfg.tapsmofrq = 30; %use 10 for 0-30Hz TFR and 30 for 30-300Hz TFR
cfg.toi = toy;
cfg.keeptrials = 'yes';
cfg.output = 'pow';

HPCspwCBD = ft_freqanalysis(cfg, ft_cbdspw);

%visualize 
cfg = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.title = 'CBD HPC below pyramidal layer sharp wave 30-300Hz';
ft_colormap(hot)
ft_singleplotTFR(cfg, HPCspwCBD)

%% to do left:
%% finding co-occurance of sharp waves and ripples
%check venn diagram and extended venn diagram to find co occurances
%% categorize into 1 2 3

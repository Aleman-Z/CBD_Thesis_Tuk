clear
clc
% %Open Ephys
% addpath(uigetdir('C:\','Select Open Ephys functions folder'));

%Fieldtrip
addpath(genpath('C:\Users\students\Documents\Tugdual\GitHub\fieldtrip'));

%CorticoHippocampal
addpath(genpath('C:\Users\students\Documents\Tugdual\GitHub\CorticoHippocampal'));

%ADRITOOLS
addpath(genpath('C:\Users\students\Documents\Tugdual\GitHub\ADRITOOLS'));


% Dowmsampling 'reference channels' and accelerometer data

acq_fhz = 30e3 ;
ds_fhz = 600;

% Which rat
opts.Resize = 'on';
opts.WindowStyle = 'modal';
opts.Interpreter = 'tex';
prompt=strcat('\bf Select a rat#. Options:','{ }',['2 ' '3 ' '4 ' '5 ' '9 ' '10 ' '11 ' '12 ' '13 ' '14 ' '15 ' '16']);
answer = inputdlg(prompt,'Input',[2 30],{''},opts);
Rat=str2double(answer{1});
          
% Select file to downsample
[file,path]=uigetfile('C:\','Select data file to downsample')
patID = strcat([path,file]);

% Which Brain Area is it?
BA = inputdlg({'Brain area recorded, or type ACC if accelerometer data'},...
              'Type your selection', [1 30]);
          
% Load
[Data, ~, ~] = load_open_ephys_data_faster(patID);

% Downsampling
 Wn = [ds_fhz/acq_fhz ]; % Cutoff=fs_new/2 Hz. 
 [b,a] = butter(3,Wn); %Filter coefficients for LPF.
 Data_filt=filtfilt(b,a,Data);
 Data_filt=downsample(Data_filt,acq_fhz/ds_fhz);
 
 patID=uigetdir('C:\','Select folder where to save')
 
 if contains (BA, 'PFC') 
     save([patID '\PFC_' file '.mat']) 
 end

 if contains (BA, 'HPC') 
     save([patID '\HPC_' file '.mat']) 
 end
 
 if contains (BA, 'ACC')
     save([patID '\ACC_' file '.mat'])
 end
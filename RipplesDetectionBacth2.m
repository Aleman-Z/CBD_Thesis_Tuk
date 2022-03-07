% RIPPLES DETECTION AND ANALYSIS
% New batch - april 2020
clear all
clc

%Fieldtrip
addpath('C:\Users\students\Documents\Tugdual\GitHub\fieldtrip');

%CorticoHippocampal
addpath('C:\Users\students\Documents\Tugdual\GitHub\CorticoHippocampal');
addpath('C:\Users\students\Documents\Tugdual\GitHub\CorticoHippocampal\Ripple_detection');

%CBD
addpath('C:\Users\students\Documents\Tugdual\GitHub\CBD');

%ADRITOOLS
addpath('C:\Users\students\Documents\Tugdual\GitHub\ADRITOOLS');

addpath('C:\Users\students\Documents\Tugdual\GitHub\bluewhitered');

addpath('C:\Users\students\Documents\Tugdual\GitHub\analysis-tools');
%% Load and downsample data

%Rat 204 - CBD

patIDhpc = 'F:\rat\cannabis\acutes batch 2\204\2020-04-03_11-24-46\102_CH43.continuous';
patIDpfc = 'F:\rat\cannabis\acutes batch 2\204\2020-04-03_11-24-46\102_CH19.continuous';

acq_fhz = 30e3;
ds_fhz = 600;
 
Wn = [ds_fhz/acq_fhz ]; % Cutoff=fs_new/2 Hz. 
[b,a] = butter(3,Wn); %Filter coefficients for LPF.

[Data, ~, ~] = load_open_ephys_data_faster(patIDhpc);
[Data2, ~, ~] = load_open_ephys_data_faster(patIDpfc);
Data = filtfilt(b,a,Data);
Data2 = filtfilt(b,a,Data2);
HPC = downsample(Data,acq_fhz/ds_fhz);
PFC = downsample(Data2,acq_fhz/ds_fhz);

%% Filter

fs = ds_fhz;
e_t = 1; %s
e_samples = e_t*(fs); 
nc = floor(length(HPC)/e_samples); %Number of epochs

[filtered_NC_ref outliers_ref rem_ep_ref] = artifacts(HPC, fs, e_t, 0, 99.5);
[filtered_NC_ref2 outliers_ref2 rem_ep_ref2] = artifacts(PFC, fs, e_t, 0, 99.5); 
HPCtmp = HPC;
PFCtmp = PFC;
index = 1;
index2 = 1;
for kk=1:nc
   if ismember(kk, rem_ep_ref)
        HPCtmp(index:index+e_samples,1) = zeros(e_samples+1,1);
        index = index + e_samples;
    else
        index = index + e_samples;
   end
   if ismember(kk, rem_ep_ref2)
        PFCtmp(max(0,index2-0.5*e_samples):(index2+1.5*e_samples),1) = ones(2*e_samples+1,1) * mean(PFCtmp); %We remove an extra part (half an epoch) before and after the artifact
        index2 = index2 + e_samples;
    else
        index2 = index2 + e_samples;
   end 
end

HPC2 = HPCtmp;
PFC2 = PFCtmp;
%% Ripples detection HPC

%Band pass
Wn2 = [100/300 299/300]; % Cutoff=fs_new/2 Hz. 
[b2,a2] = butter(3,Wn2); %Filter coefficients for LPF.
HPC_filt=filtfilt(b2,a2,HPC2);
yourtimevector = (1:length(HPC_filt))/600;

tr = 20;
[S, E, M] = findRipplesLisa(HPC_filt', yourtimevector, tr, (tr)*(1/2), [] );
M_dur = NaT(1, length(M)) - NaT(1);
for kk=1:length(M)
    M_dur(kk) = duration([0 0 M(kk)]);
end

figure()
title('Ripple detection via threshold')
tiledlayout(4,1)
t1 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 length(ripple_hpc)/600]),length(yourtimevector)),HPC_filt)
hold on
% for ii = 1:length(M)   
%         hold on
% %         patch([S(ii) E(ii) E(ii) S(ii)], [-max(HPC_filt) -max(HPC_filt) max(HPC_filt) max(HPC_filt)], 'b', 'FaceAlpha', 0.5, 'LineStyle', 'none')
%         patch([duration([0 0 S(ii)]) duration([0 0 E(ii)]) duration([0 0 E(ii)]) duration([0 0 S(ii)])], [-max(HPC_filt) -max(HPC_filt) max(HPC_filt) max(HPC_filt)], 'b', 'FaceAlpha', 0.5, 'LineStyle', 'none')
% 
% end
stem(M_dur,600*ones(size(M_dur)))
t2 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 length(ripple_hpc)/600]),length(yourtimevector)), HPC2)
title('HPC signal')
t3 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 length(ripple_hpc)/600]),length(yourtimevector)), PFC2)
title('PFC signal')
t4 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 length(ripple_hpc)/600]),length(yourtimevector)), PFC_filt)
linkaxes([t1 t2 t3 t4], 'x')


%Characterization

ripple_hpc = zeros(1,length(HPC_filt));

index = 1;
index2 = 1;
while index < length(yourtimevector) & index2<length(E)
    if S(index2) == yourtimevector(index)
        while E(index2) ~= yourtimevector(index)
            ripple_hpc(index) = 1;
            index = index + 1;
        end
        index2 = index2 + 1;
    else
        index = index + 1;
    end
end

ripple_count = sum(data2ep(ripple_hpc, 60*15, 600));
states = load('C:\Users\students\Documents\Tugdual\PCA_Scorer_batch_2\Rat204\states'); %1 rem, 2nrem, 3int
states = states.clusters;
NREM = 10*double(states == 2); %one data point = 10 sec
NREM_INT = 10*double(states ~= 1); 
tmp_nrem = sum(reshape(NREM(1:2610), 6*15,29));
tmp_nrem_int = sum(reshape(NREM_INT(1:2610), 6*15,29));

cum_ripple_hpc = cumsum(ripple_hpc);

figure
plot(linspace(duration([0 0 0]),duration([0 0 length(ripple_hpc)/600]),length(HPC_filt)), cum_ripple_hpc/600)
xlabel('Time')
ylabel('Ripple time(s)')

hold on

plot(linspace(duration([0 0 0]),duration([0 0 length(ripple_hpc)/600]),length(ripple_count)), ripple_count/600)
legend('Cumulative ripple time', 'Ripple time per 15min bin');


%normalized ripple duration over NREM duration

figure
plot(linspace(duration([0 0 0]),duration([0 0 length(ripple_hpc)/600]),length(ripple_count)), (ripple_count/600)./(tmp_nrem))
hold on
plot(linspace(duration([0 0 0]),duration([0 0 length(ripple_hpc)/600]),length(ripple_count)), (ripple_count/600)./(tmp_nrem_int))
title('Normalized ripple duration per 15min bin')
xlabel('Time')
legend('Normalized with NREM time', 'Normalizedwith NREM and intermediate time')

ripple_dur = E - S;

figure
histogram(ripple_dur)
xlabel('Duration (s)')
title('Distribution of the duration of ripples')

%% Ripples detection PFC

%Band pass
PFC_filt = filtfilt(b2,a2,PFC2);
yourtimevector = (1:length(PFC_filt))/600;
%Filter out the big spikes that are obviously artifacts (all located at close to removed artifacts)
Data_filt3 = Data_filt2;
Data_filt3(abs(Data_filt3)>=17) = 0;
Data_filt3 = Data_filt3*10;
tr2=120;
[S2, E2, M2] = findRipplesLisa2020(Data_filt3', yourtimevector, tr2, (tr2)*(1/2), [] );

figure()
plot(yourtimevector,Data_filt3)
hold on
for ii = 1:length(M2)   
        hold on
        patch([S2(ii) E2(ii) E2(ii) S2(ii)], [-max(Data_filt3) -max(Data_filt3) max(Data_filt3) max(Data_filt3)], 'b', 'FaceAlpha', 0.5, 'LineStyle', 'none')
    end
stem(M2,300*ones(size(M2)))

%% Rat 203 - vehicle


patIDhpc2031 = 'F:\rat\cannabis\acutes batch 2\203\2020-04-02_12-21-50\102_CH59.continuous'; %pyramidal layer
patIDhpc2032 = 'F:\rat\cannabis\acutes batch 2\203\2020-04-02_12-21-50\102_CH56.continuous'; %above pyramidal layer
patIDhpc2033 = 'F:\rat\cannabis\acutes batch 2\203\2020-04-02_12-21-50\102_CH48.continuous'; %below pyramidal layer
patIDpfc2031 = 'F:\rat\cannabis\acutes batch 2\203\2020-04-02_12-21-50\102_CH7.continuous'; %channel 5
patIDpfc2032 = 'F:\rat\cannabis\acutes batch 2\203\2020-04-02_12-21-50\102_CH19.continuous'; %channel 5 from end

acq_fhz = 30e3;
ds_fhz = 600;
 
Wn = [ds_fhz/acq_fhz ]; % Cutoff=fs_new/2 Hz. 
[b,a] = butter(3,Wn); %Filter coefficients for LPF.

HPC203 = [];
PFC203 = [];
pathsHPC203 = {patIDhpc2031,patIDhpc2032,patIDhpc2033};
for ii=1:3
    [Data, ~, ~] = load_open_ephys_data_faster(pathsHPC203{ii});
    Data = filtfilt(b,a,Data);
    HPC203(:,ii) = downsample(Data,acq_fhz/ds_fhz);
end

pathsPFC203 = {patIDpfc2031,patIDpfc2032};
for ii=1:2
    [Data, ~, ~] = load_open_ephys_data_faster(pathsPFC203{ii});
    Data = filtfilt(b,a,Data);
    PFC203(:,ii) = downsample(Data,acq_fhz/ds_fhz);
end

TOT = abs(HPC203(:,2)) + abs(HPC203(:,1)) + abs(HPC203(:,3)) + abs(PFC203(:,1))+ abs(PFC203(:,2));

L = length(HPC203(:,1));

figure
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), TOT)

tr = 5670; %Visual threshold
outliers = false(L,1);
index = 1;
while index<L
    if TOT(index)>=tr
        outliers((index-300):index+1999) = ones(2300,1);
        index = index+2000;
    else
        index = index +1;
    end
end

%Filter out artifacts

HPC203(outliers,:) = mean(median(HPC203));
PFC203(outliers,:) = mean(median(PFC203));


figure
tiledlayout(5,1)
tt1 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), HPC203(:,2))
tt2 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), HPC203(:,1))
tt3 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), HPC203(:,3))
tt4 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), PFC203(:,1))
tt5 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), PFC203(:,2))

linkaxes([tt1 tt2 tt3 tt4 tt5], 'x', 'y')

%% Power spectrum in the pyramidal layer
    figure
for ii=1:3
    epdata = data2ep(HPC203(:,ii), 3, 600);
    [pxx,f] = pmtm(epdata,4,[],600);
    pxx_tot{ii} = pxx;
    px=mean(pxx,2);
    s = semilogy(f,(px).*f,'LineWidth',2);
    hold on
end
    grid on
    title('Power spectrum')
    xlabel('Frequency (Hz)')
    ylabel('Spectral Power')
    legend('Pyramidal layer','Above','Below')
%% Sharp wave ripples detection

%Take two hours from minute 290

HPC203cut = HPC203(290*60*600:410*60*600,:);
PFC203cut = PFC203(290*60*600:410*60*600,:);
Lcut = length(HPC203cut(:,1));

%Band pass filter to better see the sharp waves
[c,d] = butter(3, [1/300 20/300]);
filtblwpyr = filtfilt(c,d,HPC203cut(:,3));
spw = double(filtblwpyr <= -4*std(filtblwpyr));
dspw = abs(diff(spw));

% % Detect sharp waves in the HPC layer below pyramidal layer (BPL)
% spw = double(HPC203cut(:,3) <= -600);
% dspw = abs(diff(spw));

%Get the local peak of each sharp wave
spwpeak = nan(Lcut,1);
index = 1;
while index < Lcut
    if dspw(index) == 1
        index2 = index;
        index = index+1;
        while dspw(index) == 0
            index = index + 1;
        end
        locmin = min(HPC203cut((index2+1):(index+1),3));
        indmin = find(HPC203cut(:,3) == locmin);
        spwpeak(indmin) = locmin;
        index = index+1;
    else
        index = index+1;
    end
end

% Sharp Waves characterization
N_SW = sum(spwpeak); %Number of SW
sw_epoch = ~isnan(spwpeak);
sw_epoch = reshape(sw_epoch(1:floor(Lcut/(60*600))*60*600), 60*600,floor(Lcut/(60*600))); %separate in 1min bins
sw_fhz = sum(sw_epoch)/60; %Nb of SW per second for each 1min bin
figure
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, Lcut/600, 0, 'Format', 'hh:mm:ss.SSS'),120), sw_fhz)
title('Sharp Wave rate for every 1min bin')


%Ripple detection 

[e,f] = butter(3, [90/300 200/300]);
ripplespyr = filtfilt(e,f,HPC203cut(:,1));

yourtimevector = (1:length(ripplespyr))/600;

thresh = mean(ripplespyr) + 4*std(ripplespyr);
[S, E, M] = findRipplesLisa(ripplespyr', yourtimevector, thresh, (thresh)*(1/2), [] );

M_dur = NaT(1, length(M)) - NaT(1);
for kk=1:length(M)
    M_dur(kk) = duration([0 0 M(kk)]);
end

ripple_hpc203 = zeros(1,length(HPC203cut(:,3)));

index = 1;
index2 = 1;
while index < length(yourtimevector) & index2<length(E)
    if S(index2) == yourtimevector(index)
        while E(index2) ~= yourtimevector(index)
            ripple_hpc203(index) = 1;
            index = index + 1;
        end
        index2 = index2 + 1;
    else
        index = index + 1;
    end
end

% Plots

figure
tiledlayout(4,1)
t1 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, Lcut/600, 0, 'Format', 'hh:mm:ss.SSS'),Lcut), HPC203cut(:,2))
t2 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, Lcut/600, 0, 'Format', 'hh:mm:ss.SSS'),Lcut), HPC203cut(:,1))
t3 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, Lcut/600, 0, 'Format', 'hh:mm:ss.SSS'),Lcut), ripplespyr)
hold on
stem(M_dur,600*ones(size(M_dur)))
t4 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, Lcut/600, 0, 'Format', 'hh:mm:ss.SSS'),Lcut), HPC203cut(:,3))
hold on
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, Lcut/600, 0, 'Format', 'hh:mm:ss.SSS'),Lcut), filtblwpyr)
% hold on
scatter(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, Lcut/600, 0, 'Format', 'hh:mm:ss.SSS'),Lcut),spwpeak, 'red', 'filled')

linkaxes([t1 t2 t3 t4], 'x')

%% SWR PCA

%The idea is to cut the signal in 100 (or 200) ms bins, then compute
%features for each bin. NB : 100ms = 167 values

wl = ceil(1/600 * 1000 * 200);
tmp = floor(length(HPC203cut(:,1))/wl);

bin_PL = reshape(HPC203cut(1:tmp*wl,1), wl, tmp); %Pyramidal Layer
bin_BPL = reshape(HPC203cut(1:tmp*wl,3), wl, tmp); % Below Pyramidal Layer
bin_ripples = reshape(ripple_hpc203(1:tmp*wl), wl, tmp); 
bin_filt_rip = reshape(ripplespyr(1:tmp*wl), wl, tmp); % Filtered pyramidal layer for ripples
tmpsw = ~isnan(spwpeak);
bin_sw = reshape(tmpsw(1:tmp*wl), wl, tmp);

sumbin_sw = sum(bin_sw);
sumbin_rip = sum(bin_ripples);

venn = zeros(1,tmp);
for kk=1:tmp
    venn(kk) = 1*(sumbin_sw(kk)>0 & sumbin_rip(kk) == 0) + 2*(sumbin_sw(kk) == 0 & sumbin_rip(kk) > 0) + 3*(sumbin_sw(kk) > 0 & sumbin_rip(kk) > 0);
end

C = categorical(venn,[1 2 3],{'SW without ripple','Ripple without SW', 'SWR'});
figure
histogram(C,'BarWidth',0.5)

%We only run the PCA on bins that have at least a SW or a ripple (so as not
%to have a lot of random epochs

bin_PL2 = bin_PL(:,venn~=0);
bin_BPL2 = bin_BPL(:,venn~=0);
bin_ripples2 = bin_ripples(:,venn~=0); 
bin_filt_rip2 = bin_filt_rip(:,venn~=0);
tmp2 = sum(venn~=0);
[pxx,f] = pmtm(bin_PL2,4,[],600);


featuresPCA = [];

for kk=1:tmp2
    featuresPCA(1,kk) = peak2peak(bin_BPL2(:,kk)); %For sharp waves, the larger the sharp wave, the larger the value
    featuresPCA(2, kk) = log(spectral_power(90, 200, kk, pxx, f)); %Ripple power
    featuresPCA(3, kk) = sum(bin_ripples2(:,kk)); %Ripple duration
    featuresPCA(4,kk) = max(bin_filt_rip2(:,kk)); %Amplitude of the ripple
    featuresPCA(5, kk) = log(spectral_power(0.5, 1.5, kk, pxx, f)); 
    featuresPCA(6, kk) = log(spectral_power(2, 4, kk, pxx, f)); 
    featuresPCA(7, kk) = log(spectral_power(4, 8, kk, pxx, f)); 
    featuresPCA(8, kk) = log(spectral_power(8, 20, kk, pxx, f)); 
    featuresPCA(9, kk) = log(spectral_power(30, 50, kk, pxx, f)); 
    featuresPCA(10, kk) = log(spectral_power(50, 80, kk, pxx, f)); 
    featuresPCA(11,kk) = meanfreq(pxx(:,kk), 600, [90 200]);% Mean frequency of the ripple
end


featuresPCA(:,featuresPCA(2,:)<-24) = [];
featuresPCA_norm = featuresPCA ./ max(featuresPCA')';

[coeff,~,latent,~,explained] = pca(featuresPCA_norm');
PC1 = sum(featuresPCA_norm .* coeff(:,1), 1);
PC2 = sum(featuresPCA_norm .* coeff(:,2), 1);   

[clusters, C, sumD] = kmeans([PC1' PC2'], 4);

clus_tot = zeros(1,Lcut);
clus_tot(venn~=0) = clusters;
clus_tot = clus_tot+1;
figure
% scatter(featuresPCA_norm(1,:), featuresPCA_norm(2,:))
% hold on
scatter(PC1, PC2, [], clusters)

tic
    figure
    plot(HPC203cut(:,3))
    N = wl;
    colors = ['w','b','r','g', 'y'];
    index1 = 1;
    index2 = 1;
%     tot_clusters = [];
    while index2<Lcut
        if clus_tot(index2)~= clus_tot(index2+1) || (index2 + 1) == Lcut
%                 tot_clusters((1+(index1-1)*N):index2*N) = clusters(index2);
                cl = clus_tot(index2);
                patch([(1+(index1-1)*N) index2*N index2*N (1+(index1-1)*N)], [-max(abs(HPC203cut(:,3))) -max(abs(HPC203cut(:,3))) max(abs(HPC203cut(:,3))) max(abs(HPC203cut(:,3)))], colors(cl), 'FaceAlpha', 0.2, 'LineStyle', 'none')
                hold on
                index1 = index2 +1;
        end
            index2 = index2 + 1;

    end
    
%     STATES{ii} = tot_clusters;
    title('Clustered ripples ')
    ylabel('Amplitude')
    xlabel('Time')
    toc
    
%% SWR PCA

%Other test, but rather than bins of 200ms, we gather by event

wl = ceil(1/600 * 1000 * 200);
tmp = floor(length(HPC203cut(:,1))/wl);

bin_PL = reshape(HPC203cut(1:tmp*wl,1), wl, tmp); %Pyramidal Layer
bin_BPL = reshape(HPC203cut(1:tmp*wl,3), wl, tmp); % Below Pyramidal Layer
bin_ripples = reshape(ripple_hpc203(1:tmp*wl), wl, tmp); 
bin_filt_rip = reshape(ripplespyr(1:tmp*wl), wl, tmp); % Filtered pyramidal layer for ripples
tmpsw = ~isnan(spwpeak);
bin_sw = reshape(tmpsw(1:tmp*wl), wl, tmp);

sumbin_sw = sum(bin_sw);
sumbin_rip = sum(bin_ripples);

venn = zeros(1,tmp);
for kk=1:tmp
    venn(kk) = 1*(sumbin_sw(kk)>0 & sumbin_rip(kk) == 0) + 2*(sumbin_sw(kk) == 0 & sumbin_rip(kk) > 0) + 3*(sumbin_sw(kk) > 0 & sumbin_rip(kk) > 0);
end

C = categorical(venn,[1 2 3],{'SW without ripple','Ripple without SW', 'SWR'});
figure
histogram(C,'BarWidth',0.5)

%Now we gather together the non-zero successive elements in venn
event_index = []; %Wil contain the index of beginning/end of each event
index = 1;
while index<=tmp
    if venn(index)~=0
        start = index;
        index = index+1;
        while venn(index)~=0 & index <=tmp
            index = index+1;
        end
        event_index = [event_index [start; index-1]];
    else
        index = index+1;
    end
end


%We only run the PCA every "event"

featuresPCA = zeros(5,length(event_index));

for kk=1:length(event_index)
    featuresPCA(1,kk) =  peak2peak(HPC203cut((1+wl*event_index(1,kk)-wl):(wl*event_index(2,kk)),3)); %For sharp waves, the larger the sharp wave, the larger the value
    [pxx,f] = pmtm(HPC203cut((1+wl*event_index(1,kk)-wl):(wl*event_index(2,kk)),3),4,[],600);
    featuresPCA(2, kk) = log(spectral_power(90, 200, 1, pxx, f)); %Ripple power
    featuresPCA(3, kk) = sum(sum(bin_ripples(:,event_index(1,kk):event_index(2,kk)))); %Ripple duration
    featuresPCA(4,kk) = max(max(bin_filt_rip(:,event_index(1,kk):event_index(2,kk)))); %Amplitude of the ripple
%     featuresPCA(5, kk) = log(spectral_power(0.5, 1.5, 1, pxx, f)); 
%     featuresPCA(6, kk) = log(spectral_power(2, 4, 1, pxx, f)); 
%     featuresPCA(7, kk) = log(spectral_power(4, 8, 1, pxx, f)); 
%     featuresPCA(8, kk) = log(spectral_power(8, 20, 1, pxx, f)); 
%     featuresPCA(9, kk) = log(spectral_power(30, 50, 1, pxx, f)); 
%     featuresPCA(10, kk) = log(spectral_power(50, 80, 1, pxx, f)); 
%     featuresPCA(11,kk) = double(sum(ripple_hpc203((1+wl*event_index(1,kk)-wl):(wl*event_index(2,kk))))>0)*meanfreq(pxx, 600, [90 200]);% Mean frequency of the ripple in the event (if any, 0 otherwise)
    featuresPCA(5,kk) = sum(tmpsw((1+wl*event_index(1,kk)-wl):(wl*event_index(2,kk)))); %Number of sharp waves in the event
end


featuresPCA_norm = featuresPCA ./ max(featuresPCA')';

[coeff,~,latent,~,explained] = pca(featuresPCA_norm');
PC1 = sum(featuresPCA_norm .* coeff(:,1), 1);
PC2 = sum(featuresPCA_norm .* coeff(:,2), 1);   

PC1 = PC1 ./ max(abs(PC1));
PC2 = PC2 ./ max(abs(PC2));

[clusters, C, sumD] = kmeans([PC1' PC2'], 4);

clus_tot = zeros(1,Lcut);
for kk=1:length(event_index)
    clus_tot((1+wl*event_index(1,kk)-wl):(wl*event_index(2,kk))) = clusters(kk);
end
clus_tot = clus_tot+1;
figure
% scatter(featuresPCA_norm(1,:), featuresPCA_norm(2,:))
% hold on
colors = zeros(length(clusters),3);
colors(clusters == 1,:) = repmat([0 0 1], length(colors(clusters == 1)),1);
colors(clusters == 2,:) = repmat([1 0 0], length(colors(clusters == 2)),1);
colors(clusters == 3,:) = repmat([0 1 0], length(colors(clusters == 3)),1);
colors(clusters == 4,:) = repmat([0 0 0], length(colors(clusters == 4)),1);
% colors(clusters == 5,:) = repmat([1 0 1], length(colors(clusters == 5)),1);
scatter(PC1, PC2, [], colors)

tic
    figure
    tiledlayout(2,1)

    t1 = nexttile;
    plot(HPC203cut(:,3))
    N = wl;
    colors = ['w','b','r','g', 'k', 'm'];
    index1 = 1;
    index2 = 1;
    while index2<Lcut
        if clus_tot(index2)~= clus_tot(index2+1) || (index2 + 1) == Lcut
                cl = clus_tot(index2);
                patch([index1 index2 index2 index1], [-max(abs(HPC203cut(:,3))) -max(abs(HPC203cut(:,3))) max(abs(HPC203cut(:,3))) max(abs(HPC203cut(:,3)))], colors(cl), 'FaceAlpha', 0.2, 'LineStyle', 'none')
                hold on
                index1 = index2 +1;
        end
            index2 = index2 + 1;

    end
    
    title('Clustered ripples ')
    ylabel('Amplitude')
    xlabel('Time')
    
    t2 = nexttile;
    plot(ripplespyr)
    index1 = 1;
    index2 = 1;
    while index2<Lcut
        if clus_tot(index2)~= clus_tot(index2+1) || (index2 + 1) == Lcut
                cl = clus_tot(index2);
                patch([index1 index2 index2 index1], [-max(abs(ripplespyr)) -max(abs(ripplespyr)) max(abs(ripplespyr)) max(abs(ripplespyr))], colors(cl), 'FaceAlpha', 0.2, 'LineStyle', 'none')
                hold on
                index1 = index2 +1;
        end
            index2 = index2 + 1;

    end
    ylabel('Amplitude')
    xlabel('Time')
    xlim([0 length(ripplespyr)])

    linkaxes([t1 t2], 'x')
toc


% Count of each cluster
count_clusters = [];
c1 = double(clus_tot==2);
c2 = double(clus_tot==3);
c3 = double(clus_tot==4);
c4 = double(clus_tot==5);
count_clusters(:,1) = sum(reshape(c1(1:floor(Lcut/(60*600))*60*600), 60*600,floor(Lcut/(60*600)))); %separate in 1min bins
count_clusters(:,2) = sum(reshape(c2(1:floor(Lcut/(60*600))*60*600), 60*600,floor(Lcut/(60*600))));
count_clusters(:,3) = sum(reshape(c3(1:floor(Lcut/(60*600))*60*600), 60*600,floor(Lcut/(60*600))));
count_clusters(:,4) = sum(reshape(c4(1:floor(Lcut/(60*600))*60*600), 60*600,floor(Lcut/(60*600))));

figure
plot(count_clusters(:,1)/wl, 'b')
hold on
plot(count_clusters(:,2)/wl, 'r')
hold on
plot(count_clusters(:,3)/wl, 'g')
hold on
plot(count_clusters(:,4)/wl, 'k')
title('Count of each cluster')

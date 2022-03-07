% SWR types

%% ||||| SWRS TYPES |||||

%Fieldtrip
addpath('C:\Users\students\Documents\Tugdual\GitHub\fieldtrip');

%CorticoHippocampal
addpath('C:\Users\students\Documents\Tugdual\GitHub\CorticoHippocampal');
addpath('C:\Users\students\Documents\Tugdual\GitHub\CorticoHippocampal\Ripple_detection');

%CBD
addpath('C:\Users\students\Documents\Tugdual\GitHub\CBD');

%ADRITOOLS
addpath('C:\Users\students\Documents\Tugdual\GitHub\ADRITOOLS');

 
% Additionnal package
addpath('C:\Users\students\Documents\Tugdual\GitHub\bluewhitered'); %Actually I'm not using this one here, it's for a blue/red colormap
addpath('C:\Users\students\Documents\Tugdual\GitHub\analysis-tools');


%% Rat 204 - vehicle

patIDhpc2041 = 'E:\rat\cannabis\acutes batch 2\204\2020-04-03_11-24-46\102_CH43.continuous'; %pyramidal layer
patIDhpc2042 = 'E:\rat\cannabis\acutes batch 2\204\2020-04-03_11-24-46\102_CH56.continuous'; %above pyramidal layer
patIDhpc2043 = 'E:\rat\cannabis\acutes batch 2\204\2020-04-03_11-24-46\102_CH61.continuous'; %below pyramidal layer
patIDpfc2041 = 'E:\rat\cannabis\acutes batch 2\204\2020-04-03_11-24-46\102_CH7.continuous'; %channel 5 in depth
patIDpfc2042 = 'E:\rat\cannabis\acutes batch 2\204\2020-04-03_11-24-46\102_CH19.continuous'; %channel 5 from end

acq_fhz = 30e3;
ds_fhz = 600;
 
Wn = [ds_fhz/acq_fhz ]; % Cutoff=fs_new/2 Hz. 
[b,a] = butter(3,Wn); %Filter coefficients for LPF.

HPC204 = [];
PFC204 = [];
pathsHPC204 = {patIDhpc2041,patIDhpc2042,patIDhpc2043};
for ii=1:3
    [Data, ~, ~] = load_open_ephys_data_faster(pathsHPC204{ii});
    Data = filtfilt(b,a,Data);
    HPC204(:,ii) = downsample(Data,acq_fhz/ds_fhz);
end

pathsPFC204 = {patIDpfc2041,patIDpfc2042};
for ii=1:2
    [Data, ~, ~] = load_open_ephys_data_faster(pathsPFC204{ii});
    Data = filtfilt(b,a,Data);
    PFC204(:,ii) = downsample(Data,acq_fhz/ds_fhz);
end

TOT = abs(HPC204(:,2)) + abs(HPC204(:,1)) + abs(HPC204(:,3)) + abs(PFC204(:,1))+ abs(PFC204(:,2));

L = length(HPC204(:,1));

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

HPC204(outliers,:) = mean(median(HPC204));
PFC204(outliers,:) = mean(median(PFC204));


figure
tiledlayout(5,1)
tt1 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), HPC204(:,2))
title('HPC - above pyramidal layer')
tt2 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), HPC204(:,1))
title('HPC - pyramidal layer')
tt3 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), HPC204(:,3))
title('HPC - below pyramidal layer')
tt4 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), PFC204(:,1))
title('PFC - shallow layer')
tt5 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), PFC204(:,2))
title('PFC - deep layer')

linkaxes([tt1 tt2 tt3 tt4 tt5], 'x', 'y')

%% Sharp wave ripples detection

L = length(HPC204(:,1));

%Band pass filter to better see the sharp waves
[c,d] = butter(3, [1/300 20/300]);
filtblwpyr = filtfilt(c,d,HPC204(:,3));
spw = double(filtblwpyr <= -5*std(filtblwpyr));
dspw = abs(diff(spw));

% % Detect sharp waves in the HPC layer below pyramidal layer (BPL)

%Get the local peak of each sharp wave
spwpeak = nan(L,1);
index = 1;
while index < L
    if dspw(index) == 1
        index2 = index;
        index = index+1;
        while dspw(index) == 0
            index = index + 1;
        end
        locmin = min(HPC204((index2+1):(index+1),3));
        indmin = find(HPC204(:,3) == locmin);
        spwpeak(indmin) = locmin;
        index = index+1;
    else
        index = index+1;
    end
end

% Sharp Waves characterization
N_SW = sum(spwpeak); %Number of SW
sw_epoch = ~isnan(spwpeak);
sw_epoch = reshape(sw_epoch(1:floor(L/(60*600))*60*600), 60*600,floor(L/(60*600))); %separate in 1min bins
sw_fhz = sum(sw_epoch); %Nb of SW per second for each 1min bin
figure
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),length(sw_fhz)), sw_fhz)
title('Sharp Wave rate for every 1min bin')


%Ripple detection 

[e,f] = butter(3, [90/300 200/300]);
ripplespyr = filtfilt(e,f,HPC204(:,1));

yourtimevector = (1:length(ripplespyr))/600;

thresh = mean(ripplespyr) + 5*std(ripplespyr);
[S, E, M] = findRipplesLisa(ripplespyr', yourtimevector, thresh, (thresh)*(1/2), [] );

M_dur = NaT(1, length(M)) - NaT(1);
for kk=1:length(M)
    M_dur(kk) = duration([0 0 M(kk)]);
end

ripple_hpc204 = zeros(1,length(HPC204(:,3)));

index = 1;
index2 = 1;
while index < length(yourtimevector) & index2<length(E)
    if S(index2) == yourtimevector(index)
        while E(index2) ~= yourtimevector(index)
            ripple_hpc204(index) = 1;
            index = index + 1;
        end
        index2 = index2 + 1;
    else
        index = index + 1;
    end
end

%For the Venn diagram, group by bins of 200ms
wl = 600*0.2; %200ms bins ==> 600 points = 1sec
tmp = floor(length(HPC204(:,1))/wl);

bin_PL = reshape(HPC204(1:tmp*wl,1), wl, tmp); %Pyramidal Layer
bin_BPL = reshape(HPC204(1:tmp*wl,3), wl, tmp); % Below Pyramidal Layer
bin_ripples = reshape(ripple_hpc204(1:tmp*wl), wl, tmp); 
bin_filt_rip = reshape(ripplespyr(1:tmp*wl), wl, tmp); % Filtered pyramidal layer for ripples
tmpsw = ~isnan(spwpeak);
bin_sw = reshape(tmpsw(1:tmp*wl), wl, tmp);

sumbin_sw = sum(bin_sw);
sumbin_rip = sum(bin_ripples);

venn = zeros(1,tmp);
for kk=1:tmp
    venn(kk) = 1*(sumbin_sw(kk)>0 & sumbin_rip(kk) == 0) + 2*(sumbin_sw(kk) == 0 & sumbin_rip(kk) > 0) + 3*(sumbin_sw(kk) > 0 & sumbin_rip(kk) > 0);
end

COUNT = zeros(1,3);

% Need to adjust venn for consecutive "events" = one longer event
index1 = 1;
while index1 < length(venn)
    if venn(index1) ~=0 % Beginning of an "event"
       index2 = index1+1;
       while venn(index2) ~=0 & index2 < length(venn) %We look the length of that event (or consecutive events)
           index2 = index2+1;
       end
       event = venn(index1:(index2 -1));
       if length(event) >= 2 %If that event last longer than 2 200ms bins, then it potentially is a mix of different "types" of events
           % Several cases
           if ismember(3, event)
               venn(index1:(index2-1)) = 3; %If one 200 ms bin is a SWR, then the whole event is a SWR
           elseif ismember(1, event) & ismember(2, event) 
                venn(index1:(index2-1)) = 3; %If it contains both a RIP without SW and a SW without RIP, then it's a SWR
           end
       end
       COUNT = COUNT + ismember(1,venn(index1:(index2-1)))*[1 0 0] + ismember(2,venn(index1:(index2-1)))*[0 1 0]+ ismember(3,venn(index1:(index2-1)))*[0 0 1];
       index1 = index2;
    else
        index1 = index1 +1;
    end
end


figure
b = bar([0.5 2 3.5], COUNT);
b.FaceColor = 'flat';
b.EdgeColor = 'none';
xticks([0.5 2 3.5])
xticklabels({'SW without ripple','Ripple without SW', 'SWR'})
xtickangle(25)
b.CData(1,:) = [0    0.4471    0.7412];
b.CData(2,:) = [0.9216    0.2863    0.2863];
b.CData(3,:) = [0.4667    0.7804    0.3529];
title('Number of each event type for rat 204')

% Plots

figure
tiledlayout(4,1)
% t1 = nexttile;
colorsvenn = ['b', 'r', 'g'];
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), HPC204(:,1))
% hold on
% for ii=1:length(venn)
%     if venn(ii) ~=0
%         patch([(1+(ii-1)*wl) (ii*wl) (ii*wl) (1+(ii-1)*wl)]*1/(600*24*3600), [-max(abs(HPC204(:,1))) -max(abs(HPC204(:,1))) max(abs(HPC204(:,1))) max(abs(HPC204(:,1)))], colorsvenn(venn(ii)), 'FaceAlpha', 0.5, 'LineStyle', 'none' )
%         hold on
%     end
% end
% title('Pyramidal Layer')
t1 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), ripplespyr)
hold on
stem(M_dur,600*ones(size(M_dur)))
title('Pyramidal layer filtered for ripples')
t2 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), HPC204(:,3))
hold on
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), filtblwpyr)
% hold on
scatter(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L),spwpeak, 'red', 'filled')
title('Below pyramidal layer filtered with sharp waves')
t3 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), PFC204(:,1))
hold on
for ii=1:length(venn)
    if venn(ii) ~=0
        patch([(1+(ii-1)*wl) (ii*wl) (ii*wl) (1+(ii-1)*wl)]*1/(600*24*3600), [-max(abs(HPC204(:,1))) -max(abs(HPC204(:,1))) max(abs(HPC204(:,1))) max(abs(HPC204(:,1)))], colorsvenn(venn(ii)), 'FaceAlpha', 0.5, 'LineStyle', 'none' )
        hold on
    end
end
title('PFC channel')
t4 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), [0; diff(PFC204(:,1))])
linkaxes([t1 t2 t3 t4], 'x')

% Amount of each event type vs time
events_time = [];
bins = 30/0.2; %30 seconds / 0.2sec (one datapoint is 200ms)
for ii =1:3
    tmp = venn==ii;
    tmp2 = reshape(tmp(1:floor(length(tmp)/bins)*bins), bins, floor(length(tmp)/bins));
    events_time = [events_time; sum(tmp2)];
end

figure
colors_events = [0    0.4471    0.7412; 0.9216    0.2863    0.2863; 0.4667    0.7804    0.3529];
for ii = 1:3
plot(events_time(ii,:), 'Color', colors_events(ii,:))
hold on
end

%% Amplitudes of sharp waves vs ripples
% filtblwpyr for the sharp waves
% ripplespyr for the ripples
amplitudes = [];
for ii=1:length(venn)
   if venn(ii)~=0 
       sw_vec = filtblwpyr(((ii-1)*120+1):ii*120);
       rip_vec = ripplespyr(((ii-1)*120+1):ii*120);
       pfc1 = PFC204(((ii-1)*120+1):ii*120,1);
       pfc2 = PFC204(((ii-1)*120+1):ii*120,2);
       ampl_sw = mean(filtblwpyr)-min(sw_vec);
       ampl_rip = peak2peak(rip_vec);
       ampl_pfc1 = peak2peak(pfc1);
       ampl_pfc2 = peak2peak(pfc2);
       amplitudes = [amplitudes; [ampl_sw ampl_rip ampl_pfc1 ampl_pfc2]];
   end
end

figure
scatter(amplitudes(:,1), amplitudes(:,2))
title('Amplitude of sharp waves vs ripples')
xlabel('Sharp wave amplitude') 
ylabel('Ripple amplitude')

%% REMINDER
% Venn: 1=SW, blue, 2=rip, red, 3=SWR, green
%% Spectral power in PFC after events 

epdata = data2ep(PFC204(:,1), 0.2, 600); %200ms epochs
[pxx,f] = pmtm(epdata,4,[],ds_fhz);
    
%Heatmap of power, for visualizing purposes
heatmap_data = [];
for jj=1:size(epdata,2)
   fhz = [0:2.5:200];
   for kk = 1:(length(fhz)-1)
                heatmap_data(jj,kk) = spectral_power(fhz(kk), fhz(kk+1), jj, pxx, f);
    end
end

%Correct for the artifacts that create weird powers (because signal equal
%to 0)
out = log(mean(heatmap_data'))<-2;
heatmap_data(out,:) = repelem(min(heatmap_data(~out,:))',1,sum(out))';

figure
tiledlayout(2,1)
t1 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),length(heatmap_data)), mean(heatmap_data'))
    ylim([0 200])
    hold on
    timerange = [duration([0 0 0]) duration([0 0 L/600])];
    imagesc(datenum(timerange),[min(fhz) max(fhz)], log(heatmap_data'))
    colormap(t1, 'parula')
    %colorbar();
    set(gca,'YDir','normal')
    title('Spectrogram')
    
t2 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), PFC204(:,1))
hold on
for ii=1:length(venn)
    if venn(ii) ~=0
        patch([(1+(ii-1)*wl) (ii*wl) (ii*wl) (1+(ii-1)*wl)]*1/(600*24*3600), [-max(abs(PFC204(:,1))) -max(abs(PFC204(:,1))) max(abs(PFC204(:,1))) max(abs(PFC204(:,1)))], colorsvenn(venn(ii)), 'FaceAlpha', 0.5, 'LineStyle', 'none' )
        hold on
    end
end
title('PFC channel superficial')

linkaxes([t1 t2], 'x')

%% Average over a 1sec (?) window for every event type

window_size = 10;

av_pfc_power_rip = zeros(window_size,80);
av_pfc_power_sw = zeros(window_size,80);
av_pfc_power_swr = zeros(window_size,80);
n_sw = 0;
n_rip = 0;
n_swr = 0;

index1 = 1;
while index1 < length(venn)
    if venn(index1) ~=0 % Beginning of an "event"
       index2 = index1+1;
       while venn(index2) ~=0 & index2 < length(venn) %We look at the length of that event (or consecutive events)
           index2 = index2+1;
       end
       event = venn(index1:(index2 -1));
       window_low = index1 - floor((window_size - length(event))/2);
       window_high = index2 + ceil((window_size - length(event))/2) - 1;
       pfc_spectrogram_window = heatmap_data(window_low:window_high,:);
       if event(1) == 1;
           av_pfc_power_sw = av_pfc_power_sw + pfc_spectrogram_window;
           n_sw = n_sw + 1;
       elseif event(1) == 2;
           av_pfc_power_rip = av_pfc_power_rip + pfc_spectrogram_window;
           n_rip = n_rip + 1;
       elseif event(1) == 3;
           av_pfc_power_swr = av_pfc_power_swr + pfc_spectrogram_window;
           n_swr = n_swr + 1;
       end
       index1 = index2;
    else
        index1 = index1 +1;
    end
end
av_pfc_power_sw = av_pfc_power_sw / n_sw;
av_pfc_power_rip = av_pfc_power_rip / n_rip;
av_pfc_power_swr = av_pfc_power_swr / n_swr;

figure
tiledlayout(1,3)
t1 = nexttile();
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, window_size*0.2, 0, 'Format', 'hh:mm:ss.SSS'),window_size), mean(av_pfc_power_rip'))
    ylim([0 200])
    hold on
    timerange = [duration([0 0 0]) duration([0 0 window_size*0.2])];
    imagesc(datenum(timerange),[min(fhz) max(fhz)], log(av_pfc_power_rip'))
    colormap(t1, 'parula')
    set(gca,'YDir','normal')
    title('Averaged spectrogram of PFC activity around ripples')
t2 = nexttile();
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, 0.2*window_size, 0, 'Format', 'hh:mm:ss.SSS'),window_size), mean(av_pfc_power_sw'))
    ylim([0 200])
    hold on
    timerange = [duration([0 0 0]) duration([0 0 window_size*0.2])];
    imagesc(datenum(timerange),[min(fhz) max(fhz)], log(av_pfc_power_sw'))
    colormap(t1, 'parula')
    set(gca,'YDir','normal')
    title('Averaged spectrogram of PFC activity around SWs')
t3 = nexttile();
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, 0.2*window_size, 0, 'Format', 'hh:mm:ss.SSS'),window_size), mean(av_pfc_power_swr'))
    ylim([0 200])
    hold on
    timerange = [duration([0 0 0]) duration([0 0 window_size*0.2])];
    imagesc(datenum(timerange),[min(fhz) max(fhz)], log(av_pfc_power_swr'))
    colormap(t1, 'parula')
    %colorbar();
    set(gca,'YDir','normal')
    title('Averaged spectrogram of PFC activity around SWRs')
    
    linkaxes([t1 t2 t3], 'y')
    
%% Average over a 1sec window for every event type - split low/high fhz

%0-20 Hz
figure
tiledlayout(1,3)
t1 = nexttile();
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, window_size*0.2, 0, 'Format', 'hh:mm:ss.SSS'),window_size), mean(av_pfc_power_rip'))
    ylim([0 17.5])
    hold on
    timerange = [duration([0 0 0]) duration([0 0 window_size*0.2])];
    imagesc(datenum(timerange),[fhz(1) fhz(8)], log(av_pfc_power_rip(:,1:8)'))
    colormap(t1, 'parula')
    set(gca,'YDir','normal')
    title('Averaged spectrogram of PFC activity around ripples')
t2 = nexttile();
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, 0.2*window_size, 0, 'Format', 'hh:mm:ss.SSS'),window_size), mean(av_pfc_power_sw'))
    ylim([0 17.5])
    hold on
    timerange = [duration([0 0 0]) duration([0 0 window_size*0.2])];
    imagesc(datenum(timerange),[fhz(1) fhz(8)], log(av_pfc_power_sw(:,1:8)'))
    colormap(t1, 'parula')
    set(gca,'YDir','normal')
    title('Averaged spectrogram of PFC activity around SWs')
t3 = nexttile();
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, 0.92*window_size, 0, 'Format', 'hh:mm:ss.SSS'),window_size), mean(av_pfc_power_swr'))
    ylim([0 17.5])
    hold on
    timerange = [duration([0 0 0]) duration([0 0 window_size*0.2])];
    imagesc(datenum(timerange),[fhz(1) fhz(8)], log(av_pfc_power_swr(:,1:8)'))
    colormap(t1, 'parula')
    %colorbar();
    set(gca,'YDir','normal')
    title('Averaged spectrogram of PFC activity around SWRs')
    
    linkaxes([t1 t2 t3], 'y')
    
    
% 25 - 200 Hz

figure
tiledlayout(1,3)
t1 = nexttile();
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, window_size*0.2, 0, 'Format', 'hh:mm:ss.SSS'),window_size), mean(av_pfc_power_rip'))
    ylim([25 200])
    hold on
    timerange = [duration([0 0 0]) duration([0 0 window_size*0.2])];
    imagesc(datenum(timerange),[fhz(11) fhz(80)], log(av_pfc_power_rip(:,11:80)'))
    colormap(t1, 'parula')
    set(gca,'YDir','normal')
    title('Averaged spectrogram of PFC activity around ripples')
t2 = nexttile();
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, 0.2*window_size, 0, 'Format', 'hh:mm:ss.SSS'),window_size), mean(av_pfc_power_sw'))
    ylim([25 200])
    hold on
    timerange = [duration([0 0 0]) duration([0 0 window_size*0.2])];
    imagesc(datenum(timerange),[fhz(11) fhz(80)], log(av_pfc_power_sw(:,11:80)'))
    colormap(t1, 'parula')
    set(gca,'YDir','normal')
    title('Averaged spectrogram of PFC activity around SWs')
t3 = nexttile();
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, 0.2*window_size, 0, 'Format', 'hh:mm:ss.SSS'),window_size), mean(av_pfc_power_swr'))
    ylim([25 200])
    hold on
    timerange = [duration([0 0 0]) duration([0 0 window_size*0.2])];
    imagesc(datenum(timerange),[fhz(11) fhz(80)], log(av_pfc_power_swr(:,11:80)'))
    colormap(t1, 'parula')
    %colorbar();
    set(gca,'YDir','normal')
    title('Averaged spectrogram of PFC activity around SWRs')
    
    linkaxes([t1 t2 t3], 'y')
    
    
%% DIfference between powers


%0-200 Hz
figure
tiledlayout(1,3)
t1 = nexttile();
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, window_size*0.2, 0, 'Format', 'hh:mm:ss.SSS'),window_size), mean(av_pfc_power_rip'))
    ylim([0 200])
    hold on
    timerange = [duration([0 0 0]) duration([0 0 window_size*0.2])];
    imagesc(datenum(timerange),[min(fhz) max(fhz)], log(abs(av_pfc_power_rip'- av_pfc_power_sw')))
    colormap(t1, 'parula')
    set(gca,'YDir','normal')
    title('Averaged spectrogram of PFC activity - Difference rip/sw')
t2 = nexttile();
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, 0.2*window_size, 0, 'Format', 'hh:mm:ss.SSS'),window_size), mean(av_pfc_power_sw'))
    ylim([0 200])
    hold on
    timerange = [duration([0 0 0]) duration([0 0 window_size*0.2])];
    imagesc(datenum(timerange),[min(fhz) max(fhz)], log(abs(av_pfc_power_rip'- av_pfc_power_swr')))
    colormap(t1, 'parula')
    set(gca,'YDir','normal')
    title('Averaged spectrogram of PFC activity  - Difference rip/swr')
t3 = nexttile();
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, 0.92*window_size, 0, 'Format', 'hh:mm:ss.SSS'),window_size), mean(av_pfc_power_swr'))
    ylim([0 200])
    hold on
    timerange = [duration([0 0 0]) duration([0 0 window_size*0.2])];
    imagesc(datenum(timerange),[min(fhz) max(fhz)], log(abs(av_pfc_power_swr' - av_pfc_power_sw')))
    colormap(t1, 'parula')
    %colorbar();
    set(gca,'YDir','normal')
    title('Averaged spectrogram of PFC activity  - Difference sw/swr')
    
    linkaxes([t1 t2 t3], 'y')
    

%% Spectral power in PFC after events - with the other PFC channel

epdata = data2ep(PFC204(:,2), 0.2, 600); 
[pxx,f] = pmtm(epdata,4,[],ds_fhz);
    
%Heatmap of power, for visualizing purposes
heatmap_data = [];
for jj=1:size(epdata,2)
   fhz = [0:2.5:200];
   for kk = 1:(length(fhz)-1)
                heatmap_data(jj,kk) = spectral_power(fhz(kk), fhz(kk+1), jj, pxx, f);
    end
end

%Correct for the artifacts that create weird powers (because signal equal
%to 0)
out = log(mean(heatmap_data'))<-2;
heatmap_data(out,:) = repelem(min(heatmap_data(~out,:))',1,sum(out))';

figure
tiledlayout(2,1)
t1 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),length(heatmap_data)), mean(heatmap_data'))
    ylim([0 200])
    hold on
    timerange = [duration([0 0 0]) duration([0 0 L/600])];
    imagesc(datenum(timerange),[min(fhz) max(fhz)], log(heatmap_data'))
    colormap(t1, 'parula')
    %colorbar();
    set(gca,'YDir','normal')
    title('Spectrogram')
    
t2 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), PFC204(:,2))
hold on
for ii=1:length(venn)
    if venn(ii) ~=0
        patch([(1+(ii-1)*wl) (ii*wl) (ii*wl) (1+(ii-1)*wl)]*1/(600*24*3600), [-max(abs(PFC204(:,2))) -max(abs(PFC204(:,2))) max(abs(PFC204(:,2))) max(abs(PFC204(:,2)))], colorsvenn(venn(ii)), 'FaceAlpha', 0.5, 'LineStyle', 'none' )
        hold on
    end
end
title('PFC channel deep')

linkaxes([t1 t2], 'x')

%% Looking at features of sw/rip
venn_extended = [repelem(venn,120) zeros(1,L-length(venn)*120)];
ripples_charac = [];
index1 = 1
while index1 < L
    if ripple_hpc204(index1) == 1
        index2 = index1;
        while ripple_hpc204(index2) == 1 & index2 < L
            index2 = index2 + 1;
        end
        AMPL = peak2peak(ripplespyr(index1:(index2-1)));
        [pxx,f] = pmtm(ripplespyr(index1:(index2-1)),4,[],600);
        FREQ = meanfreq(pxx, 600, [70 290]);
        DUR = (index2-index1)*1000/600;
        TYPE = venn_extended(index1);
        ripples_charac = [ripples_charac; AMPL DUR FREQ TYPE];
        index1 = index2;
    else
      index1 = index1 + 1;  
    end
end

% Amplitude of ripples

figure
h1 = histogram(ripples_charac(ripples_charac(:,4) == 3, 1),'BinWidth', 10 ,'FaceAlpha', 0.4, 'FaceColor','g');
hold on
h2 = histogram(ripples_charac(ripples_charac(:,4) == 2, 1), 'BinWidth', h1.BinWidth, 'FaceAlpha', 0.6, 'FaceColor','r');
title('Distribution of the ripples amplitudes')
xlabel('Peak to trough amplitude')
legend({'SWRs', 'Ripples only'})

% Duration of ripples
figure
h1 = histogram(ripples_charac(ripples_charac(:,4) == 3, 2),'BinWidth', 10 ,'FaceAlpha', 0.4, 'FaceColor','g');
hold on
h2 = histogram(ripples_charac(ripples_charac(:,4) == 2, 2), 'BinWidth', h1.BinWidth, 'FaceAlpha', 0.6, 'FaceColor','r');
title('Distribution of the ripples duration')
xlabel('Duration (ms)')
legend({'SWRs', 'Ripples only'})

% Mean frequency of ripples
figure
h1 = histogram(ripples_charac(ripples_charac(:,4) == 3, 3),'BinWidth', 3 ,'FaceAlpha', 0.4, 'FaceColor','g');
hold on
h2 = histogram(ripples_charac(ripples_charac(:,4) == 2, 3), 'BinWidth', h1.BinWidth, 'FaceAlpha', 0.6, 'FaceColor','r');
title('Distribution of the ripples mean frequency')
xlabel('Frequency (Hz)')
legend({'SWRs', 'Ripples only'})

SW_charac = [];
for ii = 1:L
    if ~isnan(spwpeak(ii)) 
        AMPL = -spwpeak(ii);
        TYPE = venn_extended(ii);
        SW_charac = [SW_charac; AMPL TYPE];
    end
end
% Amplitude of sharp waves
figure
h1 = histogram(SW_charac(SW_charac(:,2) == 3, 1),'BinWidth', 50 ,'FaceAlpha', 0.4, 'FaceColor','g');
hold on
h2 = histogram(SW_charac(SW_charac(:,2) == 1, 1), 'BinWidth', h1.BinWidth, 'FaceAlpha', 0.6, 'FaceColor','b');
title('Distribution of the SW amplitudes')
xlabel('SW amplitude')
legend({'SWRs', 'SWs only'})

%% Split ripples

M_dur = NaT(1, length(M)) - NaT(1);
for kk=1:length(M)
    M_dur(kk) = duration([0 0 M(kk)]);
end

ripple_hpc204 = zeros(1,length(HPC204(:,3)));

index = 1;
index2 = 1;
while index < length(yourtimevector) & index2<length(E)
    if S(index2) == yourtimevector(index)
        while E(index2) ~= yourtimevector(index)
            ripple_hpc204(index) = 1;
            index = index + 1;
        end
        index2 = index2 + 1;
    else
        index = index + 1;
    end
end

% Slow ripples
[e,f] = butter(3, [90/300 115/300]);
rip_slow = filtfilt(e,f,HPC204(:,1));
yourtimevector = (1:length(rip_slow))/600;
thresh_slow = mean(rip_slow) + 5*std(rip_slow);
[S_slow, E_slow, M_slow] = findRipplesLisa(rip_slow', yourtimevector, thresh_slow, (thresh_slow)*(1/2), [] );

% Fast ripples
[e,f] = butter(3, [115/300 200/300]);
rip_fast = filtfilt(e,f,HPC204(:,1));
thresh_fast = mean(rip_fast) + 5*std(rip_fast);
[S_fast, E_fast, M_fast] = findRipplesLisa(rip_fast', yourtimevector, thresh_fast, (thresh_fast)*(1/2), [] );

figure
t = tiledlayout(2,1);
t1 = nexttile;
plot(rip_slow)
t2 = nexttile
plot(rip_fast)
linkaxes([t1 t2], 'x')

% Features of slow/fast ripples

ripple_slow_hpc204 = zeros(1,length(HPC204(:,3)));

index = 1;
index2 = 1;
while index < length(yourtimevector) & index2<length(E_slow)
    if S_slow(index2) == yourtimevector(index)
        while E_slow(index2) ~= yourtimevector(index)
            ripple_slow_hpc204(index) = 1;
            index = index + 1;
        end
        index2 = index2 + 1;
    else
        index = index + 1;
    end
end

ripple_fast_hpc204 = zeros(1,length(HPC204(:,3)));

index = 1;
index2 = 1;
while index < length(yourtimevector) & index2<length(E_fast)
    if S_fast(index2) == yourtimevector(index)
        while E_fast(index2) ~= yourtimevector(index)
            ripple_fast_hpc204(index) = 1;
            index = index + 1;
        end
        index2 = index2 + 1;
    else
        index = index + 1;
    end
end

ripples_slow_charac = [];
index1 = 1
while index1 < L
    if ripple_slow_hpc204(index1) == 1
        index2 = index1;
        while ripple_slow_hpc204(index2) == 1 & index2 < L
            index2 = index2 + 1;
        end
        AMPL = peak2peak(rip_slow(index1:(index2-1)));
%         [pxx,f] = pmtm(rip_slow(index1:(index2-1)),4,[],600);
%         FREQ = meanfreq(pxx, 600, [70 290]);
        DUR = (index2-index1)*1000/600;
        ripples_slow_charac = [ripples_slow_charac; AMPL DUR];
        index1 = index2;
    else
      index1 = index1 + 1;  
    end
end


ripples_fast_charac = [];
index1 = 1
while index1 < L
    if ripple_fast_hpc204(index1) == 1
        index2 = index1;
        while ripple_fast_hpc204(index2) == 1 & index2 < L
            index2 = index2 + 1;
        end
        AMPL = peak2peak(rip_fast(index1:(index2-1)));
%         [pxx,f] = pmtm(rip_fast(index1:(index2-1)),4,[],600);
%         FREQ = meanfreq(pxx, 600, [70 290]);
        DUR = (index2-index1)*1000/600;
        ripples_fast_charac = [ripples_fast_charac; AMPL DUR];
        index1 = index2;
    else
      index1 = index1 + 1;  
    end
end

power_slow_fast = [];
index1 = 1
while index1 < L
    if ripple_hpc204(index1) == 1
        index2 = index1;
        while ripple_hpc204(index2) == 1 & index2 < L
            index2 = index2 + 1;
        end
        [pxx,f] = pmtm(ripplespyr(index1:(index2-1)),4,[],600);
        slopo = spectral_power(80, 110, 1, pxx, f);
        fapo = spectral_power(110, 200, 1, pxx, f);

        power_slow_fast = [power_slow_fast; slopo fapo];
        index1 = index2;
    else
      index1 = index1 + 1;  
    end
end

%plots


% Amplitude of ripples

figure
h1 = histogram(ripples_slow_charac(:, 1),'BinWidth', 10 ,'FaceAlpha', 0.4, 'FaceColor','g');
hold on
h2 = histogram(ripples_fast_charac(:, 1), 'BinWidth', h1.BinWidth, 'FaceAlpha', 0.6, 'FaceColor','r');
title('Distribution of the ripples amplitudes')
xlabel('Peak to trough amplitude')
legend({'Slow ripples', 'Fast ripples'})

% Duration of ripples
figure
h1 = histogram(ripples_slow_charac(:, 2),'BinWidth', 10 ,'FaceAlpha', 0.4, 'FaceColor','g');
hold on
h2 = histogram(ripples_fast_charac(:, 2), 'BinWidth', h1.BinWidth, 'FaceAlpha', 0.6, 'FaceColor','r');
title('Distribution of the ripples duration')
xlabel('Duration (ms)')
legend({'Slow ripples', 'Fast ripples'})

% slow power vs fast power

figure
scatter(log(power_slow_fast(:,1)), log(power_slow_fast(:,2)))
xlabel('Log-transformed slow power (90-115Hz)')
ylabel('Log-transformed fast power (115-200Hz)')
title('Slow vs fast ripple power')

% Adrian's gaussmix_slow_fast model
% all ripples
freq1 = ripples_charac(:,3);
[cutoff1, sd1] = gaussmix_slow_fast(freq1);

%Rip alone
freq2 = ripples_charac(ripples_charac(:,4) == 2,3);
[cutoff2, sd2] = gaussmix_slow_fast(freq2);

%SWRs only
freq3 = ripples_charac(ripples_charac(:,4) == 3,3);
[cutoff3, sd3] = gaussmix_slow_fast(freq3);

% ==> Gaussian mixture model does not fit: it works but the cutoff is
% obviously wrong... then should we select it by hand?
%% Look at events during UP/DOWN
[c,d] = butter(3, [0.1/300 2/300]);
up_down = filtfilt(c,d,PFC204(:,1));
data = load('C:\Users\students\Documents\Tugdual\PCA_Scorer_batch_2\Rat204\states2');
states204 = data.clusters;

%Use 20-60Hz to detect DOWN states
epdata = data2ep(PFC204(:,1), 0.2, 600); 
[pxx,f] = pmtm(epdata,4,[],ds_fhz);

power20_60 = [];
for jj=1:size(epdata,2)
      val = spectral_power(25, 180, jj, pxx, f);
      power20_60 = [power20_60; val];
end
ampl_power20_60 = [repelem(power20_60, 120) ;zeros(L-120*length(power20_60),1)];

figure
tiledlayout(4,1)
t1 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), PFC204(:,1))
t2 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), up_down)
t3 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), ampl_power20_60)


%VISUAL THRESHOLD FOR DOWN STATES : POWER > 500

UP_DOWN_log = 2*(ampl_power20_60<500) + 1 *(ampl_power20_60>=500);
t4 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L),UP_DOWN_log)
ylim([0 3])
linkaxes([t1 t2 t3 t4], 'x');

%% Other UPDOWN detection
[c,d] = butter(3, [0.1/300 1/300]);
up_down2 = filtfilt(c,d,PFC204(:,1));
portion_nrem = up_down2(10620000:12620000);
t=1:length(portion_nrem);
x=portion_nrem;
h=hilbert(x);
phase_radians=angle(h);
de=rad2deg(phase_radians);
figure
plot(t,de)
phase_degrees=mod(de,360);
plot(t,x); hold on; plot(t,phase_degrees)

venn_cut = venn_extended(10620000:12620000);

event_phase = [];
index = 1
while index < length(venn_cut)
    if venn_cut(index) ~= 0
        index2 = index;
        while venn_cut(index2+1) == venn_cut(index2) & index2 < length(venn_cut)
            index2 = index2 + 1;
        end
        mean_phase = mean(phase_degrees(index:index2));
        event_phase = [event_phase; venn_cut(index) index2-index+1 mean_phase];
        index = index2+1;
    else
        index = index + 1;
    end
end

% Now histogram of degrees per group

hist_swr = event_phase(event_phase(:,1)==3,3);
hist_sw = event_phase(event_phase(:,1)==1,3);
hist_rip = event_phase(event_phase(:,1)==2,3);

figure
tiledlayout(1,3)
t1 = nexttile;
polarhistogram(deg2rad(hist_swr), 18, 'FaceColor','g');
title('SWR');
t2 = nexttile;
polarhistogram(deg2rad(hist_sw), 18, 'FaceColor','b');
title('SW');
t3 = nexttile;
polarhistogram(deg2rad(hist_rip), 18, 'FaceColor','r');
title('Ripples');


%% Esther's UPDOWN detection

%Not necessary after checking, the one above already gives a result really
%close

%% Analysis during UP/DOWNs

% For now we restrict to the 2 big NREM periods where UP-DOWNS oscillations
% settle nicely. Boundaries are:
bounds = zeros(1,L);
bounds(4046000:5125000) = 1;
bounds(10620000:12620000)= 1;
bounds = bounds == 1;

PFC204_cut = PFC204(:,1);
UP_DOWN_cut = UP_DOWN_log;
PFC204_cut(~bounds) = 0;
UP_DOWN_cut(~bounds) = 0;
venn_cut = [repelem(venn,120) zeros(1,L-length(venn)*120)];
venn_cut(~bounds) = 0;
venn_cut2 = venn_cut(1:120:end);
UP_DOWN_cut2 = UP_DOWN_cut(1:120:end);

figure
tiledlayout(2,1)
t1 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), PFC204_cut)
for ii=1:length(venn_cut2)
    if venn_cut2(ii) ~=0
        patch([(1+(ii-1)*wl) (ii*wl) (ii*wl) (1+(ii-1)*wl)]*1/(600*24*3600), [-max(abs(PFC204_cut)) -max(abs(PFC204_cut)) max(abs(PFC204_cut)) max(abs(PFC204_cut))], colorsvenn(venn_cut2(ii)), 'FaceAlpha', 0.5, 'LineStyle', 'none' )
        hold on
    end
end
t2 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L),UP_DOWN_cut)
for ii=1:length(venn_cut2)
    if venn_cut2(ii) ~=0
        patch([(1+(ii-1)*wl) (ii*wl) (ii*wl) (1+(ii-1)*wl)]*1/(600*24*3600), [0 0 3 3], colorsvenn(venn_cut2(ii)), 'FaceAlpha', 0.5, 'LineStyle', 'none' )
        hold on
    end
end
ylim([0 3])
linkaxes([t1 t2], 'x');

% Events during UP or DOWN

%SWR - venn = 3
swr_up_down = UP_DOWN_cut2(venn_cut2 == 3);

figure
b = bar([0.5 2.5], [sum(swr_up_down==1) sum(swr_up_down==2)]);
b.FaceColor = 'flat';
b.EdgeColor = 'none';
xticks([0.5 2.5])
xticklabels({'SWR during UP','SWR during DOWN'})
xtickangle(25)
b.CData(1,:) = [0.4667    0.7804    0.3529];
b.CData(2,:) = [0.3412    0.4902    0.1412];
title('SWR count during UP/DOWN states')

%RIP alone - venn = 2
rip_up_down = UP_DOWN_cut2(venn_cut2 == 2);

figure
b = bar([0.5 2.5], [sum(rip_up_down==1) sum(rip_up_down==2)]);
b.FaceColor = 'flat';
b.EdgeColor = 'none';
xticks([0.5 2.5])
xticklabels({'RIP during UP','RIP during DOWN'})
xtickangle(25)
b.CData(1,:) = [0.9216    0.2863    0.2863];
b.CData(2,:) = [0.6000    0.1843    0.1843];
title('Ripple count during UP/DOWN states')

%SW alone - venn = 1
sw_up_down = UP_DOWN_cut2(venn_cut2 == 1);

figure
b = bar([0.5 2.5], [sum(sw_up_down==1) sum(sw_up_down==2)]);
b.FaceColor = 'flat';
b.EdgeColor = 'none';
xticks([0.5 2.5])
xticklabels({'SW during UP','SW during DOWN'})
xtickangle(25)
b.CData(1,:) = [0.1255    0.5922    0.9020];
b.CData(2,:) = [ 0    0.4471    0.7412];
title('SW count during UP/DOWN states')

% Look at duration of DOWN states that have an event in them
mat_duration = [];
index = 1
while index < length(UP_DOWN_cut2)
    if UP_DOWN_cut2(index) == 1
        index2 = index;
        while UP_DOWN_cut2(index2) == 1 & index2 < length(UP_DOWN_cut2)
            index2 = index2 +1;
        end
        dura = index2 - index;
        swr_test = ismember(3, venn_cut2(index:(index2-1)));
        mat_duration = [mat_duration; dura swr_test];
        index = index2; 
    else
        index = index +1;
    end
end

figure
h1 = histogram(mat_duration(mat_duration(:,2) == 1, 1)','BinWidth', 1 ,'FaceAlpha', 0.4, 'FaceColor',[0.3922    0.8314    0.0745]);
hold on
h2 = histogram(mat_duration(mat_duration(:,2) == 0, 1)', 'BinWidth', h1.BinWidth, 'FaceAlpha', 0.6, 'FaceColor',[0.6510    0.6510    0.6510]);
hold on
xline(mean(mat_duration(mat_duration(:,2) == 1, 1)),'--', 'Color',[0.3922    0.8314    0.0745], 'LineWidth', 5)
hold on
xline(mean(mat_duration(mat_duration(:,2) == 0, 1)),'--', 'Color',[0.6510    0.6510    0.6510], 'LineWidth', 5)
title('Distribution of the duration of DOWN states')
xlabel('Duration (in epoch of 120 ms)')
legend({'With SWR', 'Without SWR'})


figure
boxplot([mat_duration(mat_duration(:,2) == 1, 1)' mat_duration(mat_duration(:,2) == 0, 1)'], [ones(1, length(mat_duration(mat_duration(:,2) == 1, 1))) zeros(1, length(mat_duration(mat_duration(:,2) == 0, 1)))]);
title('Distribution of the duration of DOWN states')
xticks([1 2])
xticklabels({'Without SWR', 'With SWR'})
a = get(get(gca,'children'),'children');
set(a(3), 'Color', [0.1490    0.3216    0.0275]); 
set(a(4), 'Color', [0.3804    0.3804    0.3804]); 
set(a(5), 'Color', [0.3922    0.8314    0.0745]); 
set(a(6), 'Color', [0.6510    0.6510    0.6510]); 
hold on
yline(mean(mat_duration(mat_duration(:,2) == 1, 1)),'--', 'Color',[0.3922    0.8314    0.0745], 'LineWidth', 2)
hold on
yline(mean(mat_duration(mat_duration(:,2) == 0, 1)),'--', 'Color',[0.6510    0.6510    0.6510], 'LineWidth', 2)

% Phase diagram

[c,d] = butter(3, [0.1/300 1/300]);
periodic = filtfilt(c,d,PFC204(:,1));

figure
plot(PFC204(:,1))
hold on
plot(periodic)


%% Look at cortical fast oscillations during events
[c,d] = butter(3, [100/300 299/300]);
fast_oscillations = filtfilt(c,d,PFC204(:,1));

%Problem around the removed artifacts ==> set to zero around these regions
penumbra = 1*600; %ms converted into datapoints, region around artifacts that we remove additionally
index = 1;
while index < L
    if outliers(index) == 1
        index2 = index;
        while outliers(index2) == 1
            index2 = index2+1;
        end
        fast_oscillations((index-penumbra):(index2+penumbra-1)) = 0;
        index = index2+penumbra;
    else
        index = index+1;
    end
end


%Get envelope
env = abs(hilbert(fast_oscillations));

%Cut to only long NREM
fast_osci_cut = fast_oscillations;
fast_osci_cut(PFC204_cut==0) = 0;
env_cut = env;
env_cut(PFC204_cut==0) = 0;

%Detect above threshold
thresh = 2*std(fast_oscillations(PFC204_cut~=0));
index_events = env_cut>=thresh;

%Get start, end and mid of each event
EVENTS = [];
index = 1;
while index < L;
    if index_events(index) == 1;
        index2 = index;
        while index2 < L & index_events(index2) ==1
            index2 = index2+1;
        end
        EVENTS = [EVENTS; index floor(0.5*index + 0.5*(index+2-1)) index2-1];
        index = index2;
    else
        index = index+1;
    end
end

%Postprocessing

%First, if 2 consecutives events happen within 20ms, consider it's the same
%event
gap = 600*20/1000;
index = 1;
while index < length(EVENTS)
    if (EVENTS(index+1,1) - EVENTS(index,3))  <= gap
        EVENTS(index+1,1) = EVENTS(index,1);
        EVENTS(index+1,2) = floor((EVENTS(index+1,1) + EVENTS(index+1,3))/2);
        EVENTS(index,:) = [0 0 0];
    end
    index = index+1;
end
EVENTS((EVENTS(:,1) == 0),:) = [];

% Then we remove events that are very short (less than 20 ms)
outliers = zeros(1,length(EVENTS));
for ii = 1:length(EVENTS)
    if (EVENTS(ii,3)-EVENTS(ii,2)) <= gap
        outliers(ii) = 1;
    end
end
EVENTS(outliers==1,:) = [];

events = fast_osci_cut;
events(index == 0) = NaN;
events(index == 1) = max(fast_osci_cut);

index_events2 = zeros(L,1);
for ii = 1:length(EVENTS)
    index_events2(EVENTS(ii,1):EVENTS(ii,3)) = 1;
end

figure
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), fast_osci_cut)
hold on
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), env_cut)
hold on
yline(mean(fast_osci_cut)+2*std(fast_oscillations(PFC204_cut~=0)))
hold on
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), index_events*20)
% hold on
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), index_events2*20)


% Now we have our events, we want to look at the mean frequency of each
% event

frequency_fast_osci = [];
for ii=1:length(EVENTS)
    event = fast_osci_cut(EVENTS(ii,1):EVENTS(ii,3));
    [pxx,f] = pmtm(event,4,[],600);
    FREQ = meanfreq(pxx, 600, [100 297]);
    frequency_fast_osci = [frequency_fast_osci; FREQ];
end

%Now plot distribution of these frequencies

figure
histogram(frequency_fast_osci, 'BinWidth', 2);
title('Distribution of the mean frequency of fast oscillatory events in the PFC')

%NB ==> ONLY DURING NREM, because harder to detect events in REM, or short
%NREM
%% Multiplets 
% We want here to detect ripples that are close to one another
[M_multiplets, Mx_multiplets]=findmultiplets(M);
timestamps = [M_multiplets.singlets{1} M_multiplets.doublets{1} M_multiplets.triplets{1} M_multiplets.quatruplets{1}];
% find timestamps that occur more than once, and when they occur, that way
% we can know which groups are mixed up

multip = zeros(size(ripplespyr));
for k = 1:length(timestamps)
    multip(floor(timestamps(k))*600) = 1*(k<=1104) + 2*(k>1104 & k<=1480) + 3*(k>1480 & k <=1588) + 4 * (k>1588);
end

figure
tiledlayout(2,1)
t1 = nexttile();
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), ripplespyr)

t2 = nexttile();
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), multip)
ylim([-1 5])

linkaxes([t1 t2], "x")
%% PFC spectrogram during events


%% SWR PCA

%We gather by event


%Now we gather together the non-zero successive elements in venn
event_index = []; %Will contain the index of beginning/end of each event
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
epdata = data2ep(HPC204(:,1), 0.2, 600);
[pxx,f] = pmtm(epdata,4,[],400);


featuresPCA = zeros(12,length(event_index));

for kk=1:length(event_index)
    featuresPCA(1,kk) =  peak2peak(HPC204((1+wl*event_index(1,kk)-wl):(wl*event_index(2,kk)),3)); %For sharp waves, the larger the sharp wave, the larger the value
    [pxx,f] = pmtm(HPC204((1+wl*event_index(1,kk)-wl):(wl*event_index(2,kk)),3),4,[],600);
    featuresPCA(2, kk) = log(spectral_power(90, 200, 1, pxx, f)); %Ripple power
    featuresPCA(3, kk) = sum(sum(bin_ripples(:,event_index(1,kk):event_index(2,kk)))); %Ripple duration
    tmp = bin_filt_rip(:,event_index(1,kk):event_index(2,kk));
    featuresPCA(4,kk) = peak2peak(tmp(:)); %Amplitude of the ripple
    featuresPCA(5, kk) = log(spectral_power(0.5, 1.5, 1, pxx, f)); 
    featuresPCA(6, kk) = log(spectral_power(2, 4, 1, pxx, f)); 
    featuresPCA(7, kk) = log(spectral_power(4, 8, 1, pxx, f)); 
    featuresPCA(8, kk) = log(spectral_power(8, 20, 1, pxx, f)); 
    featuresPCA(9, kk) = log(spectral_power(30, 50, 1, pxx, f)); 
    featuresPCA(10, kk) = log(spectral_power(50, 80, 1, pxx, f)); 
    featuresPCA(11,kk) = double(sum(ripple_hpc204((1+wl*event_index(1,kk)-wl):(wl*event_index(2,kk))))>0)*meanfreq(pxx, 600, [90 200]);% Mean frequency of the ripple in the event (if any, 0 otherwise)
    featuresPCA(12,kk) = sum(tmpsw((1+wl*event_index(1,kk)-wl):(wl*event_index(2,kk)))); %Number of sharp waves in the event
end


featuresPCA_norm = featuresPCA ./ max(featuresPCA')';

[coeff,~,latent,~,explained] = pca(featuresPCA_norm');
PC1 = sum(featuresPCA_norm .* coeff(:,1), 1);
PC2 = sum(featuresPCA_norm .* coeff(:,2), 1);   

PC1 = PC1 ./ max(abs(PC1));
PC2 = PC2 ./ max(abs(PC2));

[clusters, C, sumD] = kmeans([PC1' PC2'], 3);

clus_tot = zeros(1,L);
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


    figure
    tiledlayout(3,1)

    t1 = nexttile;
    N = wl;
    colors = ['w','b','r','g', 'k', 'm'];
    index1 = 1;
    index2 = 1;
    while index2<L
        if clus_tot(index2)~= clus_tot(index2+1) || (index2 + 1) == L
                cl = clus_tot(index2);
                patch([index1 index2 index2 index1], [-max(abs(HPC204(:,3))) -max(abs(HPC204(:,3))) max(abs(HPC204(:,3))) max(abs(HPC204(:,3)))], colors(cl), 'FaceAlpha', 0.5, 'LineStyle', 'none')
                hold on
                index1 = index2 +1;
        end
            index2 = index2 + 1;

    end
    plot(HPC204(:,3))
    title('Clustered ripples ')
    ylabel('Amplitude')
    xlabel('Time')
    
    t2 = nexttile;
    index1 = 1;
    index2 = 1;
    while index2<L
        if clus_tot(index2)~= clus_tot(index2+1) || (index2 + 1) == L
                cl = clus_tot(index2);
                patch([index1 index2 index2 index1], [-max(abs(ripplespyr)) -max(abs(ripplespyr)) max(abs(ripplespyr)) max(abs(ripplespyr))], colors(cl), 'FaceAlpha', 0.2, 'LineStyle', 'none')
                hold on
                index1 = index2 +1;
        end
            index2 = index2 + 1;

    end
    plot(ripplespyr)
    ylabel('Amplitude')
    xlabel('Time')
    xlim([0 length(ripplespyr)])


% Count of each cluster
count_clusters = [];
c1 = double(clus_tot==2);
c2 = double(clus_tot==3);
c3 = double(clus_tot==4);
c4 = double(clus_tot==5);
count_clusters(:,1) = sum(reshape(c1(1:floor(L/(60*600))*60*600), 60*600,floor(L/(60*600)))); %separate in 1min bins
count_clusters(:,2) = sum(reshape(c2(1:floor(L/(60*600))*60*600), 60*600,floor(L/(60*600))));
count_clusters(:,3) = sum(reshape(c3(1:floor(L/(60*600))*60*600), 60*600,floor(L/(60*600))));
count_clusters(:,4) = sum(reshape(c4(1:floor(L/(60*600))*60*600), 60*600,floor(L/(60*600))));

tmp = repelem(count_clusters,600*60,1);
count_clusters = [tmp; zeros(length(ripplespyr) - length(tmp), 4)];

t3 = nexttile;

plot(count_clusters(:,1)/wl, 'b')
hold on
plot(count_clusters(:,2)/wl, 'r')
hold on
plot(count_clusters(:,3)/wl, 'g')
hold on
plot(count_clusters(:,4)/wl, 'k')
title('Count of each cluster')

    linkaxes([t1 t2 t3], 'x')

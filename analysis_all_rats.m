clear
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
%% Load all data

% Batch 1

rats = [2 3 4 5 9 10 11 12 13 15 16];
hpc_rats = [2 3 4 5 9 10 11 12];
pfc_rats = [2 3 4 5 9 10 11 13 15 16];
type = {'--','-','-','--','-','--','--','-','-','--','--'}; % solid line '-' vehicle, dashed line '--' CBD
colors = {};
for ii=1:length(rats)
    c = [rand rand rand];
    colors = [colors c];
end
type_hpc = type(1:8);
colors_hpc = colors(1:8);
type_pfc = [type(1:7)  type(9:11)];
colors_pfc = [colors(1:7)  colors(9:11)];
HPC = struct([]);
PFC = struct([]);
pxx_hpc = struct([]);
pxx_pfc = struct([]);



for ii=1:7;
    data = load(strcat('C:\Users\students\Documents\Tugdual\PLOTS\Rat',num2str(rats(ii)),'\filtered_data'));
    data = getfield(data, 'FINAL_DATA');
    HPC{ii} = data{1};
    PFC{ii} = data{2};
    
    data_pxx = load(strcat('C:\Users\students\Documents\Tugdual\PLOTS\Rat',num2str(rats(ii)),'\powers'));
    data_pxx = getfield(data_pxx, 'pxx_tot');
    pxx_hpc{ii} = data_pxx{1};
    pxx_pfc{ii} = data_pxx{2};
    
end

% For 12, only HPC
    data = load(strcat('C:\Users\students\Documents\Tugdual\PLOTS\Rat12\filtered_data'));
    data = getfield(data, 'FINAL_DATA');
    HPC{8} = data{1};
    
    data_pxx = load(strcat('C:\Users\students\Documents\Tugdual\PLOTS\Rat12\powers'));
    data_pxx = getfield(data_pxx, 'pxx_tot');
    pxx_hpc{8} = data_pxx{1};
    
    data_f = load(strcat('C:\Users\students\Documents\Tugdual\PLOTS\Rat12\freq_PS'));
    data_f = getfield(data_f, 'f');
    f{8} = data_f;

% For 13-16, only PFC

for ii=9:11;
    data = load(strcat('C:\Users\students\Documents\Tugdual\PLOTS\Rat',num2str(rats(ii)),'\filtered_data'));
    data = getfield(data, 'FINAL_DATA');
    PFC{ii-1} = data{1};
    
    data_pxx = load(strcat('C:\Users\students\Documents\Tugdual\PLOTS\Rat',num2str(rats(ii)),'\powers'));
    data_pxx = getfield(data_pxx, 'pxx_tot');
    pxx_pfc{ii-1} = data_pxx{1};
    
    data_f = load(strcat('C:\Users\students\Documents\Tugdual\PLOTS\Rat',num2str(rats(ii)),'\freq_PS'));
    data_f = getfield(data_f, 'f');
    f{ii} = data_f;
end

% frequencies from the power spectrum, same for all rats, so we take just
% one

data_f = load(strcat('C:\Users\students\Documents\Tugdual\PLOTS\Rat2\freq_PS'));
data_f = getfield(data_f, 'f');
f = data_f;

% Scored data
scored_rats = [2,5,9,10,11];
STATES = struct([]);
EVENTS = struct([]);
TRANSITIONS = struct([]);
for jj=1:length(scored_rats)
    scorat = load(strcat('C:\Users\students\Documents\Tugdual\CBD_ephys\Rat',num2str(scored_rats(jj)),'-states.mat'));
    STATES{jj} = scorat.states;
    EVENTS{jj} = scorat.events;
    TRANSITIONS{jj} = scorat.transitions;
end

% Batch 2

%% Durations
load('C:\Users\students\Documents\Tugdual\GitHub\CBD\CBD_5.mat') % To get the duration of the trials
%Starting times for Rat2-Rat5.
t1(1) = duration([12 0 0]);
t1(2) = duration([11 53 0]);
t1(3) = duration([12 0 0]);
t1(4) = duration([11 54 0]);

%Starting times for Rat9-Rat16.
t1(5) = duration([16 20 0]);
t1(6) = duration([15 21 0]);
t1(7) = duration([15 45 0]);
t1(8) = duration([17 18 0]);
t1(9) = duration([11 30 0]);
t1(10) = duration([11 30 0]);
t1(11) = duration([11 24 0]);
t1(12) = duration([11 35 0]);


t2=t1+seconds(durations_trials(:,1))'; %Values calculated from a previous iteration.

% We remove 5 min at the end of each recording as artifact might occur in
% the end and we remove the first 30 min as there is tissue relaxation
% causing weird things to happen
t3 = t2 - duration([0 5 0]);
t1 = t1 + duration([0 30 0]);

% We need to know to how much data points these removed minutes correspond
% We only need to do it for one rat, it'll be the same number of data
% points for the other rat
REMOVED = ([]);
k = 1;
for ii=[2 3 4 5 9 10 11 12 13]
    REMOVED{ii} = load(strcat('C:\Users\students\Documents\Tugdual\PLOTS\Rat',num2str(ii),'\removed_sec'));
end

%% Superimpose spectral powers

%HPC
figure
for ii=1:length(pxx_hpc)
    px=mean(pxx_hpc{ii},2);
    semilogy(f,(px).*f,'LineWidth',2, 'linestyle',type_hpc{ii}, 'color', colors_hpc{ii});
    hold on
end
grid on
title('Power spectrums of the HPC in all rats')
xlabel('Frequency (Hz)')
ylabel('Power')
legend('Rat2', 'Rat3', 'Rat4', 'Rat5', 'Rat9', 'Rat10', 'Rat11', 'Rat12')


%PFC
figure
for ii=1:length(pxx_pfc)
    px=mean(pxx_pfc{ii},2);
    semilogy(f,(px).*f,'LineWidth',2, 'linestyle',type_pfc{ii}, 'color', colors_pfc{ii});
    hold on
end
grid on
title('Power spectrums of the PFC in all rats')
xlabel('Frequency (Hz)')
ylabel('Power')
legend('Rat2', 'Rat3', 'Rat4', 'Rat5', 'Rat9', 'Rat10', 'Rat11', 'Rat13', 'Rat15', 'Rat16')
%% Plot HPCs

figure
for ii=1:8
    plot(HPC{ii} + (ii-1) * 2100)
    hold on
end
  
title('Filtered HPC signals')
legend('Rat2 CBD', 'Rat3 VEH', 'Rat4 VEH', 'Rat5 CBD', 'Rat9 VEH', 'Rat10 CBD', 'Rat11 CBD', 'Rat12 VEH')

%% Plots PFCs

figure
for ii=1:10
    plot(PFC{ii} + (ii-1) * 2100)
    hold on
end
  
title('Filtered PFC signals')
legend('Rat2 CBD', 'Rat3 VEH', 'Rat4 VEH', 'Rat5 CBD', 'Rat9 VEH', 'Rat10 CBD', 'Rat11 CBD', 'Rat13 VEH', 'Rat15 CBD', 'Rat16 CBD')

%% Plot all

figure
%rats that have both PFC and PHC
labels = struct([]);
for ii=1:7
    labels{ii} = plot(HPC{ii} + (ii-1) * 5000, 'color', colors{ii});
    hold on
    plot(PFC{ii} + 2000 + (ii-1) * 5000, 'color', colors{ii})
    hold on
    yline(3500 + (ii-1) * 5000, '--');
    hold on
end
%Rat 12 only HPC

labels{8} = plot(HPC{8} + 7*5000, 'color', colors{8});
hold on
yline(3500 + 7*5000, '--');

%Rats 13-16 only PFC
for ii=9:11
    labels{ii} = plot(PFC{ii-1} + 2000 + (ii-1)*5000, 'color', colors{ii});
    hold on
    yline(3500 + (ii-1) * 5000, '--');
    hold on
end
hold off
legend([labels{11} labels{10} labels{9} labels{8} labels{7} labels{6} labels{5} labels{4} labels{3} labels{2} labels{1}], {'Rat16 CBD', 'Rat15 CBD', 'Rat13 VEH', 'Rat12 VEH', 'Rat11 CBD', 'Rat10 CBD', 'Rat9 VEH', 'Rat5 CBD', 'Rat4 VEH', 'Rat3 VEH', 'Rat2 CBD'}, 'FontSize', 12)
xticks([])
title({'\fontsize{16}Filtered signals for all rats','\fontsize{10}CBD data','\fontsize{10}For each rat, bottom signal HPC, top signal PFC'})


%% HPC delta power heatmap for all rats

theta_powers = struct([]);
for ii=1:length(HPC)
    epoched_data = data2ep(HPC{ii}, 10, 600);
    [pxx,fhz] = pmtm(epoched_data,4,[],600);
    th_pwr = [];
    for kk=1:size(pxx,2)
        th_pwr(kk) = spectral_power(3.5, 6, kk, pxx, fhz);
    end
    MAX = max(th_pwr);
    theta_powers{ii} = th_pwr./MAX;
end

htmp = [];
lengths = cellfun(@(x) length(x), theta_powers);
maxlength = max(lengths);

for ii=1:8
    htmp = [htmp ;[theta_powers{ii} zeros(1, maxlength - length(theta_powers{ii}))]];
end

figure
h = heatmap(htmp, 'GridVisible', 'off', 'ColorLimit', [0 0.7], 'YDisplayLabels', {'Rat2', 'Rat3', 'Rat4', 'Rat5', 'Rat9', 'Rat10', 'Rat11', 'Rat12', });
colormap(flipud(gray(256)))
cdl = h.XDisplayLabels;
h.XDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));

%% Different frequency powers

for ii=1:8
    
    epoched_data = data2ep(HPC{ii}, 5, 600);
    [pxx,fhz] = pmtm(epoched_data,4,[],600);

    theta = [];
    delta = [];
    ver_low = [];
    low = [];
    low_g = [];
    high_g = [];
    ripples = [];
    for kk=1:size(pxx,2)
       theta(kk) = spectral_power(3.5, 6, kk, pxx, fhz);
       delta(kk) = spectral_power(6, 9, kk, pxx, fhz);
       ver_low(kk) = spectral_power(0.1, 1, kk, pxx, fhz);
       low(kk) = spectral_power(0.1, 20, kk, pxx, fhz);
       ripples(kk) = spectral_power(80, 120, kk, pxx, fhz);
       low_g(kk) = spectral_power(30, 45, kk, pxx, fhz);
       high_g(kk) = spectral_power(55, 80, kk, pxx, fhz);
    end

    pwr = [ver_low; theta; delta; low; low_g; high_g; ripples];
    mx = max(pwr,[],2);
    pwr = pwr ./ mx;

    figure
    h = heatmap(pwr, 'GridVisible', 'off', 'ColorLimit', [0 0.75], 'YDisplayLabels', {'0.1-1Hz', '3.5-6Hz', '6-9Hz', '0.1-20Hz', '30-45Hz', '55-80Hz', '80-120Hz'});
    colormap(flipud(gray(256)))
    cdl = h.XDisplayLabels;
    h.XDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));

end
%% Theta - SO clustering 

CLUS_DATA = struct([]);
Y = [];
THETA = struct([]);
SO = struct([]);

for ii=1:length(HPC)
    filt_data = HPC{ii};
    tmp = data2ep(filt_data, 8, 600);   

    [pxx_clus,f_clus] = pmtm(tmp,4,[],600);
    theta_pwr = [];
    slow_oscillations = [];
    for kk=1:size(tmp,2)
        theta_pwr(kk) = spectral_power(3.5, 6, kk, pxx_clus, f_clus);
        slow_oscillations(kk) = spectral_power(0.1, 1, kk, pxx_clus, f_clus);
    end
 
%      Data_to_cluster = [log(theta_pwr)' log(slow_oscillations)'];
    Data_to_cluster = [theta_pwr' slow_oscillations'];
    if ii == 2  %Rat3 needs to be adjusted by hand, by doing 3 clusters instead of 2 (an outlier creates an artifical cluster)
        clusters_id = kmeans(Data_to_cluster,3);
        [c,d] = histc(clusters_id,unique(clusters_id));
        c(find(c == min(c))) = [];
    else
        clusters_id = kmeans(Data_to_cluster,2);
        [c,d] = histc(clusters_id,unique(clusters_id));
    end
    CLUS_DATA{ii} = clusters_id;
    THETA{ii} = theta_pwr;
    SO{ii} = slow_oscillations;
    
    Y = [Y ; c'];
end

Y2 = Y;
for ii=1:size(Y2,1)
    theta = THETA{ii};
    so = SO{ii};
    if ~(mean(theta(CLUS_DATA{ii} == 1)) > mean(theta(CLUS_DATA{ii} == 2)) || mean(so(CLUS_DATA{ii} == 1)) < mean(so(CLUS_DATA{ii} == 2))) == 1
        tmp = Y2(ii,1);
        Y2(ii,1) = Y2(ii,2);
        Y2(ii,2) = tmp;
    end
end

Y2 = Y2 ./ sum(Y2,2);
X = categorical({'Rat2 CBD','Rat3 veh','Rat4 veh', 'Rat5 CBD', 'Rat9 veh', 'Rat10 CBD', 'Rat11 CBD', 'Rat12 veh'});
X = reordercats(X, {'Rat2 CBD', 'Rat5 CBD', 'Rat10 CBD', 'Rat11 CBD', 'Rat3 veh','Rat4 veh', 'Rat9 veh', 'Rat12 veh'});
figure
bar(X,Y2, 'stacked')

%% PCA clustering - HPC only

ACC = [];
REM_NREM = [];
durations_REM = [];
states_hpc = struct([]);
figure
REM_DISTRIB = ([]);
for ii=1:length(HPC)
    epdata = data2ep(HPC{ii}, 10, 600);
    [pxx,f] = pmtm(epdata,4,[],600);
    px=mean(pxx,2);

    %Heatmap of power
    heatmap_data = [];
    
    for jj=1:size(epdata,2)
        fhz = 0:0.5:100;
        for kk = 1:(length(fhz)-1)
            heatmap_data(jj,kk) = spectral_power(fhz(kk), fhz(kk+1), jj, pxx, f);
        end
        if rem(jj,100) == 0 
            display(jj); 
        end
    end
    
%     figure
%     img = imagesc(log(heatmap_data)');
%     colorbar();
%     set(gca,'YDir','normal')
 
    % PCA on several features
    % In order: slow oscillations, delta, theta, low beta, low gamma, high
    % gamma,
    features = [0.1, 1, 1, 3, 3, 6, 10, 20, 30, 45, 55, 80, 80, 150];
    % Compute spectral power for each
    data2cluster = [];
    for jj=1:1:size(epdata,2)
        for kk=1:(length(features)/2)
            data2cluster(jj,kk) = log(spectral_power(features(2*(kk-1)+1), features(2*kk), jj, pxx, f));
        end
    end
    data2cluster(:,length(features)/2+1) = data2cluster(:,3) ./ data2cluster(:,1);  %theta/SO
    data2cluster(:,length(features)/2+2) = max(abs(epdata)); % max amplitude of each epoch
     
    %normalization over the columns ie for each feature
    MAX = max(data2cluster);
    data2cluster_norm = data2cluster ./ MAX;
    

    [coeff,~,latent,~,explained] = pca(data2cluster_norm);
    
%     figure
%     variables = {'Slow oscillations 0.1-1Hz','Delta 1-3Hz','Theta 3-6Hz','Low beta 10-20 Hz','Low gamma 30-45 Hz','High gamma 55-80 Hz', 'Ripples 80-150 Hz', 'Theta/SO ratio', 'Epoch max amplitude'};
%     img = imagesc(coeff);
%     colormap(bluewhitered(256));
%     colorbar();
%     set(gca,'YDir','normal')
%     set(gca,'Ytick',1:9,'YTickLabel',variables)
%     set(gca,'Xtick',1:9,'XTickLabel',num2str(round(explained,2)))
%     xlabel('Principal components (1st - left to 9th - right) and percentage of variance explained by them')
%     title(strcat('Weigth of each PCA feature in each principal component - Rat ', num2str(Rat)))
    
    % We retrieve the 2 first components
    PC1 = sum(data2cluster_norm' .* coeff(:,1), 1);
    PC2 = sum(data2cluster_norm' .* coeff(:,2), 1);
    % PC3 = sum(data2cluster_norm' .* coeff(:,3), 1);
    
    
%     clusters = kmeans([log(PC1)'  real(log(PC2))'],2);
    clusters = kmeans([PC1'  PC2'],3, 'MaxIter', 500);
%     figure
%     scatter(PC1, PC2, [], clusters)
%     xlabel('1st component')
%     ylabel('2nd component')
%     title(strcat('Kmeans clustering in the PC1 - PC2 space for rat ',num2str(Rat)))

    
    %If one epoch alone is into other cluster, make it into same cluster
    %than rest ==> Seems like a bad idea ?
    index = 2;
    while index <= (length(clusters)-1)
        if clusters(index-1) == clusters(index+1) & clusters(index) ~= clusters(index-1)
            clusters(index) = clusters(index-1);
            index = index+1;
        end
        index = index+1;
    end
    
    % We set REM in cluster 1, the two other clusters are NREM and
    % intermediate
    ampl = data2cluster(:,9);
    %first find the cluster corresponding to REM and NREM
    %starting values
    rems = 1;
    nrems = 1;
    minimum = mean(ampl(clusters==1));
    maximum = mean(ampl(clusters==1));
    for jj=unique(clusters)'
      if  mean(ampl(clusters==jj)) < minimum
          rems = jj;
          minimum = mean(ampl(clusters==jj));
      end
      if mean(ampl(clusters==jj)) > maximum
          nrems = jj;
          maximum = mean(ampl(clusters==jj));
      end
    end
    if length(unique(clusters))>2
        tmp = unique(clusters);
        tmp(tmp == rems) = [];
        tmp(tmp == nrems) = [];
        int = tmp;
    end
    %now that you found the cluster corresponding to rem and nrem and intermediate, put them in
    %position 1 and 2 respectively
    clusters(clusters == rems) = 97;
    clusters(clusters == nrems) = 98;
    if length(unique(clusters))>2
        clusters(clusters == int) = 3;
    end
    clusters(clusters == 97) = 1;
    clusters(clusters == 98) = 2;         
    
    N = size(epdata,1);
    tic
 
    subplot(length(HPC), 1, ii)
    colors = ['w','b','r','g'];
    index1 = 1;
    index2 = 1;
    tot_clusters = [];
    while index2<length(clusters)
        if clusters(index2)~= clusters(index2+1) || (index2 + 1) == length(clusters)
                tot_clusters((1+(index1-1)*N):index2*N) = clusters(index2);
                cl = clusters(index2);
                patch([(1+(index1-1)*N) index2*N index2*N (1+(index1-1)*N)], [-max(abs(HPC{ii})) -max(abs(HPC{ii})) max(abs(HPC{ii})) max(abs(HPC{ii}))], colors(cl), 'FaceAlpha', 0.2, 'LineStyle', 'none')
                hold on
                index1 = index2 +1;
        end
            index2 = index2 + 1;

    end
    plot(HPC{ii})
    title(strcat('Clustered signal of the HPC of rat ', num2str(hpc_rats(ii))))
    ylabel('Amplitude')
    xlabel('Time')
    toc
    
    tmp = tot_clusters;
    tmp(tmp == 3) = 2;
    states_hpc{ii} = tmp;
    
    %Scored data
    if ismember(hpc_rats(ii),[2 5 9 10 11 12]) %Rats with manual scoring data
        load(strcat('C:\Users\students\Documents\Tugdual\CBD_ephys\Rat',num2str(hpc_rats(ii)),'-states.mat'))
        scored_data = states;
        tot_scored_data = [];
        for jj=1:length(scored_data)
            tot_scored_data((1+(jj-1)*600):jj*600) = scored_data(jj);
        end
        
        removed_sec = REMOVED{hpc_rats(ii)};
        removed_sec = removed_sec.removed;
        rem_sec2 = removed_sec(2);
        %remove first 30 min
        tot_scored_data = tot_scored_data((1+rem_sec2):length(tot_scored_data));
        %the two sets, from manual scoring and from the clustering algo, start
        %at the same time, now just adjust the end to have the same length
        tot_scored_data = tot_scored_data(1:length(tot_clusters));


        %For all rats: scored data, 3 NREM, 5 REM, clustered_data, 1 REM, 2
        %NREM, 3, if there is, intermediate
        toremove = tot_scored_data == 0;
        scored = tot_scored_data(~toremove);
        algo = tot_clusters(~toremove);
        scored(scored==3) = 2;
        scored(scored==5) = 1;
        algo(algo == 3) = 2; %for the accuracy test, we only consider 2 states, so intermediate is regarded as NREM
        correct = scored == algo;
        accuracy = sum(correct)/length(correct);
        
        CM = confusionmat(scored,algo);
        acc_rem = CM(1,1) / (CM(1,1) + CM(1,2) + CM(2,1));
        acc_nrem = CM(2,2) / (CM(2,2) + CM(1,2) + CM(2,1));
        ACC(hpc_rats(ii),:) = [accuracy acc_rem acc_nrem];
    end
    %Get duration of each REM episode
    REM_dur = [];
    ind = 1;
    while ind<length(tot_clusters)
       dur = 0;
       if tot_clusters(ind) == 1
           while tot_clusters(ind) == 1 & ind<length(tot_clusters)
               dur = dur + 1;
               ind = ind + 1;
           end
           if dur >6000 %We remove the one epoch long episodes
               REM_dur = [REM_dur dur];
           else
               dur = 0;
           end
       end
       ind = ind+1;
    end
    REM_DISTRIB{ii} = REM_dur;
    durations_REM(ii,:) = [median(REM_dur) length(REM_dur)];
    REM_NREM(ii,:) = [length(tot_clusters(tot_clusters == 1 ))/length(tot_clusters) length(tot_clusters(tot_clusters ~= 1 ))/length(tot_clusters)];
end

figure
cat = {'cbd', 'veh','veh','cbd', 'veh', 'cbd', 'cbd', 'veh'};
cbd_veh_hpc = categorical(cat);
scatter(durations_REM(:,1), durations_REM(:,2), [], cbd_veh_hpc);


figure
X = categorical({'Rat2 CBD','Rat3 veh','Rat4 veh', 'Rat5 CBD', 'Rat9 veh', 'Rat10 CBD', 'Rat11 CBD', 'Rat12 veh'});
X = reordercats(X, {'Rat2 CBD', 'Rat5 CBD', 'Rat10 CBD', 'Rat11 CBD', 'Rat3 veh','Rat4 veh', 'Rat9 veh', 'Rat12 veh'});
bar(X, REM_NREM, 'stacked')

figure
plot_acc = ACC(ACC(:,1) ~= 0, :);
Y = categorical({'Rat2 CBD', 'Rat5 CBD', 'Rat9 veh', 'Rat10 CBD', 'Rat11 CBD', 'Rat12 veh'});
X = categorical({'Global accuracy', 'REM accuracy', 'NREM accuracy'});
heatmap(plot_acc, 'XData', X, 'YData', Y);
% set(gca, 'XTick', 1:3, 'XTickLabel', X)
% set(gca, 'YTick', 1:6, 'YTickLabel', Y)
colorbar()
title('Accuracy of the PCA clustering compared to manual clustering')

figure
colors = ['r','b','b','r','b','r','r','b'];
for ii=1:8
    subplot(4,2,ii)
    histogram(REM_DISTRIB{ii}/600,20, 'FaceColor', colors(ii))
    title(strcat('Distribution of the REM bouts durations for rat', num2str(hpc_rats(ii)),'-',cat(ii)))
%     pd = fitdist(REM_DISTRIB{ii},'Normal')   
%     pdS = fitdist(REM_DISTRIB{ii}','Kernel');
%     x = linspace(0,1.5*max(REM_DISTRIB{ii}),100);
%     y = pdf(pdS,x);
%     semilogx(x,y)
    hold on
end
% legend('Rat2 CBD','Rat3 veh','Rat4 veh', 'Rat5 CBD', 'Rat9 veh', 'Rat10 CBD', 'Rat11 CBD', 'Rat12 veh')
%% PCA clustering - PFC only

ACC2 = [];
REM_NREM2 = [];
durations_REM2 = [];
states_pfc = struct([]);
figure
REM_DISTRIB2 = ([]);
for ii=1:length(PFC)
    epdata = data2ep(PFC{ii}, 10, 600);
    [pxx,f] = pmtm(epdata,4,[],600);
    px=mean(pxx,2);

    %Heatmap of power
    heatmap_data = [];
    
    for jj=1:size(epdata,2)
        fhz = 0:0.5:100;
        for kk = 1:(length(fhz)-1)
            heatmap_data(jj,kk) = spectral_power(fhz(kk), fhz(kk+1), jj, pxx, f);
        end
        if rem(jj,100) == 0 
            display(jj); 
        end
    end
    
%     figure
%     img = imagesc(log(heatmap_data)');
%     colorbar();
%     set(gca,'YDir','normal')
 
    % PCA on several features
    % In order: slow oscillations, delta, theta, low beta, low gamma, high
    % gamma,
    features = [0.1, 1, 1, 3, 3, 6, 10, 20, 30, 45, 55, 80, 80, 150];
    % Compute spectral power for each
    data2cluster = [];
    for jj=1:1:size(epdata,2)
        for kk=1:(length(features)/2)
            data2cluster(jj,kk) = log(spectral_power(features(2*(kk-1)+1), features(2*kk), jj, pxx, f));
        end
    end
    data2cluster(:,length(features)/2+1) = data2cluster(:,3) ./ data2cluster(:,1);  %theta/SO
    data2cluster(:,length(features)/2+2) = max(abs(epdata)); % max amplitude of each epoch
     
    %normalization over the columns ie for each feature
    MAX = max(data2cluster);
    data2cluster_norm = data2cluster ./ MAX;
    

    [coeff,~,latent,~,explained] = pca(data2cluster_norm);
    
%     figure
%     variables = {'Slow oscillations 0.1-1Hz','Delta 1-3Hz','Theta 3-6Hz','Low beta 10-20 Hz','Low gamma 30-45 Hz','High gamma 55-80 Hz', 'Ripples 80-150 Hz', 'Theta/SO ratio', 'Epoch max amplitude'};
%     img = imagesc(coeff);
%     colormap(bluewhitered(256));
%     colorbar();
%     set(gca,'YDir','normal')
%     set(gca,'Ytick',1:9,'YTickLabel',variables)
%     set(gca,'Xtick',1:9,'XTickLabel',num2str(round(explained,2)))
%     xlabel('Principal components (1st - left to 9th - right) and percentage of variance explained by them')
%     title(strcat('Weigth of each PCA feature in each principal component - Rat ', num2str(Rat)))
    
    % We retrieve the 2 first components
    PC1 = sum(data2cluster_norm' .* coeff(:,1), 1);
    PC2 = sum(data2cluster_norm' .* coeff(:,2), 1);
    % PC3 = sum(data2cluster_norm' .* coeff(:,3), 1);
    
    
%     clusters = kmeans([log(PC1)'  real(log(PC2))'],2);
    clusters = kmeans([PC1'  PC2'],3, 'MaxIter', 500);
%     figure
%     scatter(PC1, PC2, [], clusters)
%     xlabel('1st component')
%     ylabel('2nd component')
%     title(strcat('Kmeans clustering in the PC1 - PC2 space for rat ',num2str(Rat)))

    
    %If one epoch alone is into other cluster, make it into same cluster
    %than rest ==> Seems like a bad idea ?
    index = 2;
    while index <= (length(clusters)-1)
        if clusters(index-1) == clusters(index+1) & clusters(index) ~= clusters(index-1)
            clusters(index) = clusters(index-1);
            index = index+1;
        end
        index = index+1;
    end
    
    % We set REM in cluster 1, the two other clusters are NREM and
    % intermediate
    ampl = data2cluster(:,9);
    %first find the cluster corresponding to REM and NREM
    %starting values
    rems = 1;
    nrems = 1;
    minimum = mean(ampl(clusters==1));
    maximum = mean(ampl(clusters==1));
    for jj=unique(clusters)'
      if  mean(ampl(clusters==jj)) < minimum
          rems = jj;
          minimum = mean(ampl(clusters==jj));
      end
      if mean(ampl(clusters==jj)) > maximum
          nrems = jj;
          maximum = mean(ampl(clusters==jj));
      end
    end
    if length(unique(clusters))>2
        tmp = unique(clusters);
        tmp(tmp == rems) = [];
        tmp(tmp == nrems) = [];
        int = tmp;
    end
    %now that you found the cluster corresponding to rem and nrem and intermediate, put them in
    %position 1 and 2 respectively
    clusters(clusters == rems) = 97;
    clusters(clusters == nrems) = 98;
    if length(unique(clusters))>2
        clusters(clusters == int) = 3;
    end
    clusters(clusters == 97) = 1;
    clusters(clusters == 98) = 2;         
    
    N = size(epdata,1);
    tic
 
    subplot(length(PFC), 1, ii)
    colors = ['w','b','r','g'];
    index1 = 1;
    index2 = 1;
    tot_clusters = [];
    while index2<length(clusters)
        if clusters(index2)~= clusters(index2+1) || (index2 + 1) == length(clusters)
                tot_clusters((1+(index1-1)*N):index2*N) = clusters(index2);
                cl = clusters(index2);
                patch([(1+(index1-1)*N) index2*N index2*N (1+(index1-1)*N)], [-max(abs(PFC{ii})) -max(abs(PFC{ii})) max(abs(PFC{ii})) max(abs(PFC{ii}))], colors(cl), 'FaceAlpha', 0.2, 'LineStyle', 'none')
                hold on
                index1 = index2 +1;
        end
            index2 = index2 + 1;

    end
    plot(PFC{ii})
    title(strcat('Clustered signal of the PFC of rat ', num2str(pfc_rats(ii))))
    ylabel('Amplitude')
    xlabel('Time')
    toc
    
    tmp = tot_clusters;
    tmp(tmp == 3) = 2;
    states_pfc{ii} = tmp;
    
    %Scored data
    if ismember(pfc_rats(ii),[2 5 9 10 11 13]) %Rats with manual scoring data
        load(strcat('C:\Users\students\Documents\Tugdual\CBD_ephys\Rat',num2str(pfc_rats(ii)),'-states.mat'))
        scored_data = states;
        tot_scored_data = [];
        for jj=1:length(scored_data)
            tot_scored_data((1+(jj-1)*600):jj*600) = scored_data(jj);
        end
        
        removed_sec = REMOVED{pfc_rats(ii)};
        removed_sec = removed_sec.removed;
        rem_sec2 = removed_sec(2);
        %remove first 30 min
        tot_scored_data = tot_scored_data((1+rem_sec2):length(tot_scored_data));
        %the two sets, from manual scoring and from the clustering algo, start
        %at the same time, now just adjust the end to have the same length
        tot_scored_data = tot_scored_data(1:length(tot_clusters));


        %For all rats: scored data, 3 NREM, 5 REM, clustered_data, 1 REM, 2
        %NREM, 3, if there is, intermediate
        toremove = tot_scored_data == 0;
        scored = tot_scored_data(~toremove);
        algo = tot_clusters(~toremove);
        scored(scored==3) = 2;
        scored(scored==5) = 1;
        algo(algo == 3) = 2; %for the accuracy test, we only consider 2 states, so intermediate is regarded as NREM
        correct = scored == algo;
        accuracy = sum(correct)/length(correct);
        
        CM = confusionmat(scored,algo);
        acc_rem = CM(1,1) / (CM(1,1) + CM(1,2) + CM(2,1));
        acc_nrem = CM(2,2) / (CM(2,2) + CM(1,2) + CM(2,1));
        ACC2(pfc_rats(ii),:) = [accuracy acc_rem acc_nrem];
    end
    %Get duration of each REM episode
    REM_dur = [];
    ind = 1;
    while ind<length(tot_clusters)
       dur = 0;
       if tot_clusters(ind) == 1
           while tot_clusters(ind) == 1 & ind<length(tot_clusters)
               dur = dur + 1;
               ind = ind + 1;
           end
           if dur >6000 %We remove the one epoch long episodes
               REM_dur = [REM_dur dur];
           else
               dur = 0;
           end
       end
       ind = ind+1;
    end
    REM_DISTRIB2{ii} = REM_dur;
    durations_REM2(ii,:) = [median(REM_dur) length(REM_dur)];
    REM_NREM2(ii,:) = [length(tot_clusters(tot_clusters == 1 ))/length(tot_clusters) length(tot_clusters(tot_clusters ~= 1 ))/length(tot_clusters)];
end

figure
cat = {'cbd', 'veh','veh','cbd', 'veh', 'cbd', 'cbd', 'veh', 'cbd' ,'cbd'};
cbd_veh_pfc = categorical(cat);
scatter(durations_REM2(:,1), durations_REM2(:,2), [], cbd_veh_pfc);


figure
X = categorical({'Rat2 CBD','Rat3 veh','Rat4 veh', 'Rat5 CBD', 'Rat9 veh', 'Rat10 CBD', 'Rat11 CBD', 'Rat13 veh', 'Rat15 CBD', 'Rat16 CBD'});
X = reordercats(X, {'Rat2 CBD', 'Rat5 CBD', 'Rat10 CBD', 'Rat11 CBD','Rat15 CBD', 'Rat16 CBD', 'Rat3 veh','Rat4 veh', 'Rat9 veh', 'Rat13 veh'});
bar(X, REM_NREM2, 'stacked')

figure
plot_acc = ACC2(ACC2(:,1) ~= 0, :);
Y = categorical({'Rat2 CBD', 'Rat5 CBD', 'Rat9 veh', 'Rat10 CBD', 'Rat11 CBD', 'Rat13 veh'});
X = categorical({'Global accuracy', 'REM accuracy', 'NREM accuracy'});
heatmap(plot_acc, 'XData', X, 'YData', Y);
% set(gca, 'XTick', 1:3, 'XTickLabel', X)
% set(gca, 'YTick', 1:6, 'YTickLabel', Y)
colorbar()
title('Accuracy of the PCA clustering compared to manual clustering')

figure
colors = ['r','b','b','r','b','r','r','b', 'r', 'r'];
for ii=1:10
    subplot(5,2,ii)
    histogram(REM_DISTRIB2{ii}/600,20, 'FaceColor', colors(ii))
    title(strcat('Distribution of the REM bouts durations for rat', num2str(pfc_rats(ii)),'-',cat(ii)))
%     pd = fitdist(REM_DISTRIB2{ii},'Normal')   
%     pdS = fitdist(REM_DISTRIB2{ii}','Kernel');
%     x = linspace(0,1.5*max(REM_DISTRIB2{ii}),100);
%     y = pdf(pdS,x);
%     semilogx(x,y)
    hold on
end
% legend('Rat2 CBD','Rat3 veh','Rat4 veh', 'Rat5 CBD', 'Rat9 veh', 'Rat10 CBD', 'Rat11 CBD', 'Rat12 veh')

%% Comparison PCA HPC / PFC

ACC_pfchpc = [];

for ii=1:7 %Rats that habe both channels, ie 2 3 4 5 9 10 11
    acc_tot = sum(states_hpc{ii} == states_pfc{ii}) / length(states_hpc{ii});
    CM = confusionmat(states_hpc{ii}, states_pfc{ii});
    REM_acc = CM(1,1) / (CM(1,1) + CM(2,1) + CM(1,2));
    NREM_acc = CM(2,2) / (CM(2,2) + CM(2,1) + CM(1,2));
    ACC_pfchpc(ii,:) = [acc_tot REM_acc NREM_acc];
end
figure
Y = categorical({'Rat2 CBD', 'Rat3 veh', 'Rat4 veh', 'Rat5 CBD', 'Rat9 veh', 'Rat10 CBD', 'Rat11 CBD'});
X = categorical({'Global accuracy', 'REM accuracy', 'NREM accuracy'});
heatmap(ACC_pfchpc, 'XData', X, 'YData', Y)
title('Accuracy HPC PCA versus PFC PCA clusterings')

%% PCA clustering with both HPC and PFC features

ACC3 = [];
REM_NREM3 = [];
durations_REM3 = [];
states_hpcpfc = struct([]);
figure
REM_DISTRIB3 = ([]);
rats_both = [2 3 4 5 9 10 11];
for ii=1:length(rats_both)
    epdata = data2ep(HPC{ii}, 10, 600);
    epdata2 = data2ep(PFC{ii}, 10, 600);
    [pxx,f] = pmtm(epdata,4,[],600);
    [pxx2,f2] = pmtm(epdata2,4,[],600);

    

    % PCA on several features
    % In order: slow oscillations, delta, theta, low beta, low gamma, high
    % gamma,
    features = [0.1, 1, 1, 3, 3, 6, 10, 20, 30, 45, 55, 80, 80, 150];
    % Compute spectral power for each
    data2cluster = [];
    for jj=1:1:size(epdata,2)
        for kk=1:(length(features)/2)
            data2cluster(jj,kk) = log(spectral_power(features(2*(kk-1)+1), features(2*kk), jj, pxx, f));
            data2cluster(jj,kk +(length(features)/2)) = log(spectral_power(features(2*(kk-1)+1), features(2*kk), jj, pxx2, f2));
        end
    end
    data2cluster(:,size(data2cluster,2)+1) = data2cluster(:,3) ./ data2cluster(:,1);  %theta/SO HPC
    data2cluster(:,size(data2cluster,2)+1) = max(abs(epdata)); % max amplitude of each epoch HPC
    data2cluster(:,size(data2cluster,2)+1) = data2cluster(:,10) ./ data2cluster(:,8);  %theta/SO HPC
    data2cluster(:,size(data2cluster,2)+1) = max(abs(epdata)); % max amplitude of each epoch HPC
     
    %normalization over the columns ie for each feature
    MAX = max(data2cluster);
    data2cluster_norm = data2cluster ./ MAX;
    

    [coeff,~,latent,~,explained] = pca(data2cluster_norm);
    
%     figure
%     variables = {'Slow oscillations 0.1-1Hz','Delta 1-3Hz','Theta 3-6Hz','Low beta 10-20 Hz','Low gamma 30-45 Hz','High gamma 55-80 Hz', 'Ripples 80-150 Hz', 'Theta/SO ratio', 'Epoch max amplitude'};
%     img = imagesc(coeff);
%     colormap(bluewhitered(256));
%     colorbar();
%     set(gca,'YDir','normal')
%     set(gca,'Ytick',1:9,'YTickLabel',variables)
%     set(gca,'Xtick',1:9,'XTickLabel',num2str(round(explained,2)))
%     xlabel('Principal components (1st - left to 9th - right) and percentage of variance explained by them')
%     title(strcat('Weigth of each PCA feature in each principal component - Rat ', num2str(Rat)))
    
    % We retrieve the 2 first components
    PC1 = sum(data2cluster_norm' .* coeff(:,1), 1);
    PC2 = sum(data2cluster_norm' .* coeff(:,2), 1);
    % PC3 = sum(data2cluster_norm' .* coeff(:,3), 1);
    
    
%     clusters = kmeans([log(PC1)'  real(log(PC2))'],2);
    clusters = kmeans([PC1'  PC2'],3, 'MaxIter', 500);
%     figure
%     scatter(PC1, PC2, [], clusters)
%     xlabel('1st component')
%     ylabel('2nd component')
%     title(strcat('Kmeans clustering in the PC1 - PC2 space for rat ',num2str(Rat)))

    
    %If one epoch alone is into other cluster, make it into same cluster
    %than rest ==> Seems like a bad idea ?
    index = 2;
    while index <= (length(clusters)-1)
        if clusters(index-1) == clusters(index+1) & clusters(index) ~= clusters(index-1)
            clusters(index) = clusters(index-1);
            index = index+1;
        end
        index = index+1;
    end
    
    % We set REM in cluster 1, the two other clusters are NREM and
    % intermediate
    ampl = data2cluster(:,9);
    %first find the cluster corresponding to REM and NREM
    %starting values
    rems = 1;
    nrems = 1;
    minimum = mean(ampl(clusters==1));
    maximum = mean(ampl(clusters==1));
    for jj=unique(clusters)'
      if  mean(ampl(clusters==jj)) < minimum
          rems = jj;
          minimum = mean(ampl(clusters==jj));
      end
      if mean(ampl(clusters==jj)) > maximum
          nrems = jj;
          maximum = mean(ampl(clusters==jj));
      end
    end
    if length(unique(clusters))>2
        tmp = unique(clusters);
        tmp(tmp == rems) = [];
        tmp(tmp == nrems) = [];
        int = tmp;
    end
    %now that you found the cluster corresponding to rem and nrem and intermediate, put them in
    %position 1 and 2 respectively
    clusters(clusters == rems) = 97;
    clusters(clusters == nrems) = 98;
    if length(unique(clusters))>2
        clusters(clusters == int) = 3;
    end
    clusters(clusters == 97) = 1;
    clusters(clusters == 98) = 2;         
    
    N = size(epdata,1);
    tic
 
    subplot(7, 1, ii)
    colors = ['w','b','r','g'];
    index1 = 1;
    index2 = 1;
    tot_clusters = [];
    while index2<length(clusters)
        if clusters(index2)~= clusters(index2+1) || (index2 + 1) == length(clusters)
                tot_clusters((1+(index1-1)*N):index2*N) = clusters(index2);
                cl = clusters(index2);
                patch([(1+(index1-1)*N) index2*N index2*N (1+(index1-1)*N)], [-max(abs(PFC{ii})) -max(abs(PFC{ii})) max(abs(PFC{ii})) max(abs(PFC{ii}))], colors(cl), 'FaceAlpha', 0.2, 'LineStyle', 'none')
                hold on
                index1 = index2 +1;
        end
            index2 = index2 + 1;

    end
    plot(PFC{ii})
    title(strcat('Clustered signal of the PFC of rat ', num2str(pfc_rats(ii))))
    ylabel('Amplitude')
    xlabel('Time')
    toc
    
    tmp = tot_clusters;
    tmp(tmp == 3) = 2;
    states_hpcpfc{ii} = tmp;
    
    %Scored data
    if ismember(pfc_rats(ii),[2 5 9 10 11 13]) %Rats with manual scoring data
        load(strcat('C:\Users\students\Documents\Tugdual\CBD_ephys\Rat',num2str(pfc_rats(ii)),'-states.mat'))
        scored_data = states;
        tot_scored_data = [];
        for jj=1:length(scored_data)
            tot_scored_data((1+(jj-1)*600):jj*600) = scored_data(jj);
        end
        
        removed_sec = REMOVED{pfc_rats(ii)};
        removed_sec = removed_sec.removed;
        rem_sec2 = removed_sec(2);
        %remove first 30 min
        tot_scored_data = tot_scored_data((1+rem_sec2):length(tot_scored_data));
        %the two sets, from manual scoring and from the clustering algo, start
        %at the same time, now just adjust the end to have the same length
        tot_scored_data = tot_scored_data(1:length(tot_clusters));


        %For all rats: scored data, 3 NREM, 5 REM, clustered_data, 1 REM, 2
        %NREM, 3, if there is, intermediate
        toremove = tot_scored_data == 0;
        scored = tot_scored_data(~toremove);
        algo = tot_clusters(~toremove);
        scored(scored==3) = 2;
        scored(scored==5) = 1;
        algo(algo == 3) = 2; %for the accuracy test, we only consider 2 states, so intermediate is regarded as NREM
        correct = scored == algo;
        accuracy = sum(correct)/length(correct);
        
        CM = confusionmat(scored,algo);
        acc_rem = CM(1,1) / (CM(1,1) + CM(1,2) + CM(2,1));
        acc_nrem = CM(2,2) / (CM(2,2) + CM(1,2) + CM(2,1));
        ACC3(pfc_rats(ii),:) = [accuracy acc_rem acc_nrem];
    end
    %Get duration of each REM episode
    REM_dur = [];
    ind = 1;
    while ind<length(tot_clusters)
       dur = 0;
       if tot_clusters(ind) == 1
           while tot_clusters(ind) == 1 & ind<length(tot_clusters)
               dur = dur + 1;
               ind = ind + 1;
           end
           if dur >6000 %We remove the one epoch long episodes
               REM_dur = [REM_dur dur];
           else
               dur = 0;
           end
       end
       ind = ind+1;
    end
    REM_DISTRIB3{ii} = REM_dur;
    durations_REM3(ii,:) = [median(REM_dur) length(REM_dur)];
    REM_NREM3(ii,:) = [length(tot_clusters(tot_clusters == 1 ))/length(tot_clusters) length(tot_clusters(tot_clusters ~= 1 ))/length(tot_clusters)];
end

figure
cat = {'cbd', 'veh','veh','cbd', 'veh', 'cbd', 'cbd'};
cbd_veh_pfc = categorical(cat);
scatter(durations_REM3(:,1), durations_REM3(:,2), [], cbd_veh_pfc);


figure
X = categorical({'Rat2 CBD','Rat3 veh','Rat4 veh', 'Rat5 CBD', 'Rat9 veh', 'Rat10 CBD', 'Rat11 CBD'});
X = reordercats(X, {'Rat2 CBD', 'Rat5 CBD', 'Rat10 CBD', 'Rat11 CBD', 'Rat3 veh','Rat4 veh', 'Rat9 veh'});
bar(X, REM_NREM3, 'stacked')

figure
plot_acc = ACC3(ACC3(:,1) ~= 0, :);
Y = categorical({'Rat2 CBD','Rat5 CBD', 'Rat9 veh', 'Rat10 CBD', 'Rat11 CBD'});
X = categorical({'Global accuracy', 'REM accuracy', 'NREM accuracy'});
heatmap(plot_acc, 'XData', X, 'YData', Y);
% set(gca, 'XTick', 1:3, 'XTickLabel', X)
% set(gca, 'YTick', 1:6, 'YTickLabel', Y)
colorbar()
title('Accuracy of the PCA clustering compared to manual clustering')

figure
colors = ['r','b','b','r','b','r','r'];
for ii=1:7
    subplot(5,2,ii)
    histogram(REM_DISTRIB3{ii}/600,20, 'FaceColor', colors(ii))
    title(strcat('Distribution of the REM bouts durations for rat', num2str(pfc_rats(ii)),'-',cat(ii)))
%     pd = fitdist(REM_DISTRIB3{ii},'Normal')   
%     pdS = fitdist(REM_DISTRIB3{ii}','Kernel');
%     x = linspace(0,1.5*max(REM_DISTRIB3{ii}),100);
%     y = pdf(pdS,x);
%     semilogx(x,y)
    hold on
end
% legend('Rat2 CBD','Rat3 veh','Rat4 veh', 'Rat5 CBD', 'Rat9 veh', 'Rat10 CBD', 'Rat11 CBD', 'Rat12 veh')


%% Ripple detection

% HPC

% select rat ID in the HPC table  ... NB : ripples in rats 4 (3), 5 (4) and
% 10 (6)
ratID = 4;
data = HPC{ratID}; 

%Band pass
 Wn = [100/300 299/300]; % Cutoff=fs_new/2 Hz. 
 [b,a] = butter(3,Wn); %Filter coefficients for LPF.
 Data_filt=filtfilt(b,a,data);
yourtimevector = (1:length(Data_filt))/600;
tr=28;
[S, E, M] = findRipplesLisa(Data_filt', yourtimevector, tr, (tr)*(1/2), [] );

figure()
plot(yourtimevector,Data_filt)
hold on
for ii = 1:length(M)   
        hold on
        patch([S(ii) E(ii) E(ii) S(ii)], [-max(Data_filt) -max(Data_filt) max(Data_filt) max(Data_filt)], 'b', 'FaceAlpha', 0.5, 'LineStyle', 'none')
    end
stem(M,600*ones(size(M)))

% PFC

data2 = PFC{ratID}; %Make sure it is the same rat

%Band pass
 Wn = [100/300 299/300]; % Cutoff=fs_new/2 Hz. 
 [b,a] = butter(3,Wn); %Filter coefficients for LPF.
 Data_filt2=filtfilt(b,a,data2);
yourtimevector = (1:length(Data_filt2))/600;
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

%% Cluster with ripples

ripple_log = zeros(1,length(HPC{ratID}));

index = 1;
index2 = 1;
while index < length(yourtimevector) & index2<length(E)
    if S(index2) == yourtimevector(index)
        while E(index2) ~= yourtimevector(index)
            ripple_log(index) = 1;
            index = index + 1;
        end
        index2 = index2 + 1;
    else
        index = index + 1;
    end
end

epdata = data2ep(HPC{ratID}, 3, 600);
epripples = data2ep(ripple_log,3 ,600);

ripple_rate = sum(epripples, 1);

figure
bar(1:length(ripple_rate), ripple_rate)

%% Group the ripples that are close

% Get a logical vector that says which elements have non-zero data in them.
nonZeroElements = ripple_log;
% Define the closest regions can be.  If they are farther away than this,
% then they will be considered as separate regions.
minSeparation = 1800;
nonZeroElements = ~bwareaopen(~nonZeroElements, minSeparation);
[labeledRegions, numRegions] = bwlabel(nonZeroElements);

%Epoch the data
epripples_gr = data2ep(nonZeroElements, 3, 600);
ripple_rate_gr = sum(epripples_gr, 1);

figure
bar(1:length(ripple_rate), ripple_rate_gr)
hold on 
bar(1:length(ripple_rate), ripple_rate)

% Do it on the ripple count per epoch instead?

nonZeroElements2 = ripple_rate~=0;
minSeparation = 1;
nonZeroElements2 = ~bwareaopen(~nonZeroElements2, minSeparation);

figure
bar(1:length(ripple_rate), nonZeroElements2)


%% Ripples durations

dur_rip = E - S;
figure
histogram(dur_rip, 50)
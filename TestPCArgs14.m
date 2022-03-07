% Test PCA clustering on rgs14 data
clear all
clc
addpath('C:\Users\students\Documents\Tugdual\GitHub\analysis-tools');
addpath 'C:\Users\students\Documents\Tugdual\GitHub\CBD';
%% Load data - rat3

%21-22 10 2019

%Raw Data

%post1
post1 = struct([]);
filename = {'100_CH36_0', '100_CH21_0', '100_AUX1_0', '100_AUX2_0', '100_AUX3_0'};
%choose file from on of these
% pathid = 'E:\rat\Rat_OS_Ephys_RGS14\Rat_OS_Ephys_RGS14_rat3_357152\OS_Ephys_RGS14_Rat3_357152_SD6_OR_21-22_10_2019\2019-10-21_10-32-28_Post_Trial1\';
% pathid = 'E:\rat\Rat_OS_Ephys_RGS14\Rat_OS_Ephys_RGS14_rat3_357152\OS_Ephys_RGS14_Rat3_357152_SD6_OR_21-22_10_2019\2019-10-21_11-22-52_Post_Trial2\';
% pathid = 'E:\rat\Rat_OS_Ephys_RGS14\Rat_OS_Ephys_RGS14_rat3_357152\OS_Ephys_RGS14_Rat3_357152_SD6_OR_21-22_10_2019\2019-10-21_12-13-16_Post_Trial3\';
% pathid = 'E:\rat\Rat_OS_Ephys_RGS14\Rat_OS_Ephys_RGS14_rat3_357152\OS_Ephys_RGS14_Rat3_357152_SD6_OR_21-22_10_2019\2019-10-21_13-03-39_Post_Trial4\';
pathid = 'E:\rat\Rat_OS_Ephys_RGS14\Rat_OS_Ephys_RGS14_rat3_357152\OS_Ephys_RGS14_Rat3_357152_SD6_OR_21-22_10_2019\2019-10-21_13-54-07_Post_Trial5\';
for i=1:5 %in the order: hpc, pfc, auxiliary 1, 2 and 3
    file = strcat(pathid,filename(i),'.continuous');
    file = file{1};
    post1{i} = load_open_ephys_data(file);  
end
res_acc = sqrt(post1{3}.^2 + post1{4}.^2 + post1{5}.^2);
post1{6} = res_acc;

%Scored Data

%Choose accordingly
% post1_scored = load('C:\Users\students\Documents\Tugdual\rgs14Scored\Rat3\2019-10-21_10-32-28_Post_trial1-states');
% post1_scored = load('C:\Users\students\Documents\Tugdual\rgs14Scored\Rat3\2019-10-21_11-22-52_Post_trial2-states');
% post1_scored = load('C:\Users\students\Documents\Tugdual\rgs14Scored\Rat3\2019-10-21_12-13-16_Post_trial3-states');
% post1_scored = load('C:\Users\students\Documents\Tugdual\rgs14Scored\Rat3\2019-10-21_13-03-39_Post_trial4-states');
post1_scored = load('C:\Users\students\Documents\Tugdual\rgs14Scored\Rat3\2019-10-21_13-54-07_Post_trial5-states');
post1_states = post1_scored.states;
post1_transitions = post1_scored.transitions;

fs = 30e3;
%% PCA clustering with both HPC and PFC features

    e_t = 3;
    ephpc = data2ep(post1{1}, e_t, fs);
    eppfc = data2ep(post1{2}, e_t, fs);
    epaux = data2ep(post1{6}, e_t, fs);
    
    [pxxhpc,f] = pmtm(ephpc,4,[],fs);
    [pxxpfc,f2] = pmtm(eppfc,4,[],fs);

    %First detect wake with auxiliary channel
    %For this we retrieve the enveloppe of the signal, and check the
    %maximum
    
    var = (max(epaux)-min(epaux))/max(max(epaux)-min(epaux));
    
    outl = var>0.05;
    outl = double(outl);
    outl(outl==0) = 2;
%     test_clus = kmeans([(1:length(var))'/length(var) (max(epaux)'-min(epaux)')/max(max(epaux)-min(epaux))],2);
    
    m1 = mean(var(outl == 1));
    m2 = mean(var(outl == 2));
    %We manually put 1 = wake and 2 = sleep
    if m1 < m2
        outl(outl == 1) = 3;
        outl(outl == 2) = 1;
        outl(outl == 3) = 2;
    end
    
   % Smoothing
   %If sudden sleep within wakeful period, we assume it is still wake
   index = 1;
   while index < length(outl)
       if outl(index) == 1
           index = index +1;
       else
           count = 1;
           while outl(index) == 2 & index < length(outl);
               index = index +1;
               count = count +1;
           end
           if count <=10
               outl((index-count):(index-1)) = linspace(1,1,count);
           end
       end
   end
   
    figure
    plot(var);
    hold on
    for jj=1:length(outl)
        colors = ['w','b','r','g'];
        for kk = 1:length(unique(outl))     
            if outl(jj) == kk
                hold on
                patch([jj (jj+1) (jj+1) jj], [0 0 max(abs(var)) max(abs(var))], colors(kk), 'FaceAlpha', 0.2, 'LineStyle', 'none')
            end
        end
    end
    
    %We prepare the data to run PCA only on the sleep parts 
    ratio = size(ephpc,2)/size(epaux,2);
    for ii=1:floor(length(outl)/ratio)
        outl2(ii) = floor(mean(outl((ratio*(ii-1)+1):(ratio*ii))));
    end
    
    pxxhpc = pxxhpc(:,outl2==2);
    pxxpfc = pxxpfc(:,outl2==2);

    % PCA on several features
    % In order: slow oscillations, delta, theta, low gamma, high
    % gamma, ripples
    features = [0.5, 3, 3, 5, 5, 11, 30, 45, 55, 80, 80, 150];
    % Compute spectral power for each
    data2cluster = [];
    for jj=1:1:size(pxxhpc,2)
        for kk=1:(length(features)/2)
            data2cluster(jj,kk) = log(spectral_power(features(2*(kk-1)+1), features(2*kk), jj, pxxhpc, f));
            data2cluster(jj,kk +(length(features)/2)) = log(spectral_power(features(2*(kk-1)+1), features(2*kk), jj, pxxpfc, f));
        end
    end
    data2cluster(:,size(data2cluster,2)+1) = data2cluster(:,3) ./ data2cluster(:,2);  %theta/delta HPC
    data2cluster(:,size(data2cluster,2)+1) = max(abs(ephpc(:,outl2==2))); % max amplitude of each epoch HPC
    data2cluster(:,size(data2cluster,2)+1) = data2cluster(:,9) ./ data2cluster(:,8);  %theta/delta PFC
    data2cluster(:,size(data2cluster,2)+1) = max(abs(eppfc(:,outl2==2))); % max amplitude of each epoch PFC

     
    %normalization over the columns ie for each feature
    MAX = max(data2cluster);
    data2cluster_norm = data2cluster ./ MAX;
    

    [coeff,~,latent,~,explained] = pca(data2cluster_norm);
    

    % We retrieve the 2 first components
    PC1 = sum(data2cluster_norm' .* coeff(:,1), 1);
    PC2 = sum(data2cluster_norm' .* coeff(:,2), 1);
    % PC3 = sum(data2cluster_norm' .* coeff(:,3), 1);
    
    
%     clusters = kmeans([log(PC1)'  real(log(PC2))'],2);
    clusters = kmeans([PC1'  PC2'],2, 'MaxIter', 500);
    figure
    scatter(PC1, PC2, [], clusters)
    xlabel('1st component')
    ylabel('2nd component')
    title('Kmeans clustering in the PC1 - PC2 space for rat3 ')

    
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
    %position 3 and 5 respectively
    clusters(clusters == rems) = 97;
    clusters(clusters == nrems) = 98;
    if length(unique(clusters))>2
        clusters(clusters == int) = 4;
    end
    clusters(clusters == 97) = 5;
    clusters(clusters == 98) = 3;         
    
    states = outl2;
    states(states == 2) = clusters;
    
    figure
    plot(post1_states)
    hold on
    plot(repelem(states,e_t))
    
    N = size(ephpc,1);
%     tic
%  
%     colors = ['w','b','r','g'];
%     index1 = 1;
%     index2 = 1;
%     tot_clusters = [];
%     while index2<length(clusters)
%         if clusters(index2)~= clusters(index2+1) || (index2 + 1) == length(clusters)
%                 tot_clusters((1+(index1-1)*N):index2*N) = clusters(index2);
%                 cl = clusters(index2);
%                 patch([(1+(index1-1)*N) index2*N index2*N (1+(index1-1)*N)], [-max(abs(post1{1})) -max(abs(post1{1})) max(abs(post1{1})) max(abs(post1{1}))], colors(cl), 'FaceAlpha', 0.2, 'LineStyle', 'none')
%                 hold on
%                 index1 = index2 +1;
%         end
%             index2 = index2 + 1;
% 
%     end
%     plot(post1{1})
%     title('Clustered signal of the HPC of rat 3')
%     ylabel('Amplitude')
%     xlabel('Time')
%     toc
%     
%     tmp = tot_clusters;
%     tmp(tmp == 3) = 2;
%     states_hpcpfc{ii} = tmp;
%     
%     
%     %Scored data
%     if ismember(pfc_rats(ii),[2 5 9 10 11 13]) %Rats with manual scoring data
%         load(strcat('C:\Users\students\Documents\Tugdual\CBD_ephys\Rat',num2str(pfc_rats(ii)),'-states.mat'))
%         scored_data = states;
%         tot_scored_data = [];
%         for jj=1:length(scored_data)
%             tot_scored_data((1+(jj-1)*600):jj*600) = scored_data(jj);
%         end
%         
%         removed_sec = REMOVED{pfc_rats(ii)};
%         removed_sec = removed_sec.removed;
%         rem_sec2 = removed_sec(2);
%         %remove first 30 min
%         tot_scored_data = tot_scored_data((1+rem_sec2):length(tot_scored_data));
%         %the two sets, from manual scoring and from the clustering algo, start
%         %at the same time, now just adjust the end to have the same length
%         tot_scored_data = tot_scored_data(1:length(tot_clusters));
% 
% 
%         %For all rats: scored data, 3 NREM, 5 REM, clustered_data, 1 REM, 2
%         %NREM, 3, if there is, intermediate
%         toremove = tot_scored_data == 0;
%         scored = tot_scored_data(~toremove);
%         algo = tot_clusters(~toremove);
%         scored(scored==3) = 2;
%         scored(scored==5) = 1;
%         algo(algo == 3) = 2; %for the accuracy test, we only consider 2 states, so intermediate is regarded as NREM
%         correct = scored == algo;
%         accuracy = sum(correct)/length(correct);
%         
%         CM = confusionmat(scored,algo);
%         acc_rem = CM(1,1) / (CM(1,1) + CM(1,2) + CM(2,1));
%         acc_nrem = CM(2,2) / (CM(2,2) + CM(1,2) + CM(2,1));
%         ACC3(pfc_rats(ii),:) = [accuracy acc_rem acc_nrem];
%     end
%     %Get duration of each REM episode
%     REM_dur = [];
%     ind = 1;
%     while ind<length(tot_clusters)
%        dur = 0;
%        if tot_clusters(ind) == 1
%            while tot_clusters(ind) == 1 & ind<length(tot_clusters)
%                dur = dur + 1;
%                ind = ind + 1;
%            end
%            if dur >6000 %We remove the one epoch long episodes
%                REM_dur = [REM_dur dur];
%            else
%                dur = 0;
%            end
%        end
%        ind = ind+1;
%     end
%     REM_DISTRIB3{ii} = REM_dur;
%     durations_REM3(ii,:) = [median(REM_dur) length(REM_dur)];
%     REM_NREM3(ii,:) = [length(tot_clusters(tot_clusters == 1 ))/length(tot_clusters) length(tot_clusters(tot_clusters ~= 1 ))/length(tot_clusters)];
% 
% 
% figure
% cat = {'cbd', 'veh','veh','cbd', 'veh', 'cbd', 'cbd'};
% cbd_veh_pfc = categorical(cat);
% scatter(durations_REM3(:,1), durations_REM3(:,2), [], cbd_veh_pfc);
% 
% 
% figure
% X = categorical({'Rat2 CBD','Rat3 veh','Rat4 veh', 'Rat5 CBD', 'Rat9 veh', 'Rat10 CBD', 'Rat11 CBD'});
% X = reordercats(X, {'Rat2 CBD', 'Rat5 CBD', 'Rat10 CBD', 'Rat11 CBD', 'Rat3 veh','Rat4 veh', 'Rat9 veh'});
% bar(X, REM_NREM3, 'stacked')
% 
% figure
% plot_acc = ACC3(ACC3(:,1) ~= 0, :);
% Y = categorical({'Rat2 CBD','Rat5 CBD', 'Rat9 veh', 'Rat10 CBD', 'Rat11 CBD'});
% X = categorical({'Global accuracy', 'REM accuracy', 'NREM accuracy'});
% heatmap(plot_acc, 'XData', X, 'YData', Y);
% % set(gca, 'XTick', 1:3, 'XTickLabel', X)
% % set(gca, 'YTick', 1:6, 'YTickLabel', Y)
% colorbar()
% title('Accuracy of the PCA clustering compared to manual clustering')
% 
% figure
% colors = ['r','b','b','r','b','r','r'];
% for ii=1:7
%     subplot(5,2,ii)
%     histogram(REM_DISTRIB3{ii}/600,20, 'FaceColor', colors(ii))
%     title(strcat('Distribution of the REM bouts durations for rat', num2str(pfc_rats(ii)),'-',cat(ii)))
% %     pd = fitdist(REM_DISTRIB3{ii},'Normal')   
% %     pdS = fitdist(REM_DISTRIB3{ii}','Kernel');
% %     x = linspace(0,1.5*max(REM_DISTRIB3{ii}),100);
% %     y = pdf(pdS,x);
% %     semilogx(x,y)
%     hold on
% end
% % legend('Rat2 CBD','Rat3 veh','Rat4 veh', 'Rat5 CBD', 'Rat9 veh', 'Rat10 CBD', 'Rat11 CBD', 'Rat12 veh')



%  for l=1:length(S)k
ro=600;
l=1

    [sig,p,q,cont,sig_pq]=convert2fieldtrip({S},{E},{yourtimevector},{ripplespyr'},l,{M},{ripplespyr'},ro)
    p=p.';
    q=q.';
    sig_pq=sig_pq.';
    %%
    
%     for i=1:5
i=14
        plot(p{i})
        hold on
%     end

%%
timecell=create_timecell(ro,length(p))
%timecell=timecell;
freqrange=[30:10:300];
toy=[-1:.01:1];
 
ft_data1 = [];
ft_data1.fsample = 600;
ft_data1.trial = p(1,1:end); 
% ft_data1.time = (timecell(1,1:end));
ft_data1.time = timecell(1,1:end)

ft_data1
ft_data1.label = { 'PFC'};

% Compute spectrogram

cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'dpss';
cfg.foi = freqrange;
cfg.t_ftimwin = .1 * ones(size(cfg.foi));
cfg.tapsmofrq = 20;
cfg.toi=toy;
cfg.keeptrials = 'yes';
cfg.output         = 'pow';

freq = ft_freqanalysis(cfg, ft_data1);
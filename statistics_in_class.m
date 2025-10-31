clc;close all; clear;
addpath("C:\Users\rahma\Downloads\");
load sampleEEGdata.mat;



% Input params

chan2use= 'fcz';
min_freq = 3;
max_freq = 30;
num_frex = 20;


% Get TFR ; use morlet wavelet no. 6 

channel_idx = find(strcmpi({EEG.chanlocs.labels},chan2use));

srate = EEG.srate;

frex = logspace(log10(min_freq),log10(max_freq),num_frex);

for ev = 1:size(EEG.data,3)
    [cfs(:,:,ev),~]= morletWaveletTransform(EEG.data(channel_idx,:,ev),srate,frex,6,2);
end

power_all = abs(cfs).^2;

% Edge trimming 
time_s=dsearchn(EEG.times',-500);
time_e=dsearchn(EEG.times',1200);
eegpower= power_all(:,time_s:time_e,:);
tftimes = EEG.times(time_s:time_e);
nTimepoints = length(tftimes);

base_idx = [dsearchn(tftimes',-500),dsearchn(tftimes',-100)];

eegpower_baseline = squeeze(mean(power_all(:,base_idx(1):base_idx(2),:),2));
% log based baseline normalization 

% for trials = 1:size(EEG.data,3)
%     eegpower(:,:,trials) = 10*log10(eegpower(:,:,trials) ./ mean(eegpower_baseline,2));
% end
realmean= zeros(20,436,99);
for f = (1:20)
    for tr =(1:99)
        realmean(f,:,tr) = 10*log10(eegpower(f,:,tr)/eegpower_baseline(f,tr));
    end
    
end

% permute the time series using circshift 
real_avg_trial = squeeze(mean(realmean,3));

perm_trials = zeros(20,436,1000);

for per=1:1000              % 1000 permutations for each trial independently 
    shuffled_data = zeros(20,436,99);
    for tr = 1:99
        random_shift = randi(436);
        shuffled_data(:,:,tr)= circshift(realmean(:,:,tr),random_shift,2);
    end
    perm_trials(:,:,per) = mean(shuffled_data, 3);              % we store the 20 x 436 average trials of the shuffled data for 1000 permutations 
end

% calculate z score

z_map = (real_avg_trial- mean(perm_trials,3))./ std(perm_trials,0,3);

figure;
imagesc(tftimes,frex,z_map);

xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Z-score Map (Unthresholded)');
colorbar;

set(gca, 'YDir', 'normal', 'YScale', 'log', 'YTick', logspace(log10(min_freq), log10(max_freq), 5));
%thresholding 

threshold = norminv(0.975, 0, 1);

z_map_thresholded = z_map;
z_map_thresholded(abs(z_map)<threshold)= 0 ;

figure;
imagesc(tftimes,frex,z_map_thresholded);

xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Z-score Map (Thresholded)');
colorbar;

set(gca, 'YDir', 'normal', 'YScale', 'log', 'YTick', logspace(log10(min_freq), log10(max_freq), 5));


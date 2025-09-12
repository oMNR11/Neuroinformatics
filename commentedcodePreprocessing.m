%% -------------------------
% Load EEG dataset
% -------------------------

% Initialize EEGLAB. This is a crucial first step for using its functions.
eeglab

% Load a sample EEG dataset from the EEGLAB distribution.
% 'eeglab_data.set' is the filename.
% 'fullfile' is used to create a file path that works on different operating systems.
% 'fileparts(which('eeglab'))' finds the directory where EEGLAB is installed.
EEG = pop_loadset('eeglab_data.set', fullfile(fileparts(which('eeglab')), 'sample_data'));

% Extract the raw EEG data into a new variable for easier access.
% 'EEG.data' contains the EEG signal, with dimensions: channels x samples.
data = EEG.data;

% Get the labels for each channel from the 'EEG.chanlocs' structure.
chan_labels = {EEG.chanlocs.labels};

%% -------------------------
% Step 1: Filtering: Notch filter at 60 Hz
% -------------------------

% Calculate the normalized center frequency (wo) and bandwidth (bo) for a 60 Hz notch filter.
% 'wo' is the center frequency normalized by the Nyquist frequency (EEG.srate/2).
% 'bo' is the bandwidth, set to a fraction of the center frequency to control the filter's sharpness.
wo = 60/(EEG.srate/2); bo = wo/35;

% Design a second-order IIR (Infinite Impulse Response) notch filter.
% 'designNotchPeakIIR' is a function that creates the filter coefficients.
% 'CenterFrequency' specifies the frequency to be notched out.
% 'Bandwidth' determines the width of the frequency band to be attenuated.
% 'Response' specifies the filter type, in this case, a 'notch' filter.
% 'bn' and 'an' are the numerator and denominator coefficients of the filter.
[bn,an] = designNotchPeakIIR(CenterFrequency= wo, Bandwidth= bo, Response="notch");

% Apply the filter to the data using 'filtfilt' which performs zero-phase filtering.
% 'filtfilt(bn,an,data')'' applies the filter to the data transposed (') to filter along the correct dimension (samples).
% The result is then transposed back (') to restore the original channel x samples orientation.
data_filt = filtfilt(bn,an,data')';

%% -------------------------
% Step 2: Plot PSD before & after filtering
% -------------------------

% Create a new figure window.
figure;
% Plot the power spectral density (PSD) of the original, unfiltered data.
% 'spectopo' is an EEGLAB function for plotting spectra.
% The inputs are the data, number of data points, and sampling rate.
[~,~,~,~,~] = spectopo(data, size(data,2), EEG.srate);

% Create another figure window.
figure;
% Plot the PSD of the filtered data to visually inspect the removal of the 60 Hz noise.
[~,~,~,~,~] = spectopo(data_filt, size(data_filt,2), EEG.srate);

%% -------------------------
% Step 3: Referencing: Common Average Reference
% -------------------------

% Calculate the common average reference (CAR).
% This finds the mean of the voltage across all channels at each time point.
% 'mean(data_filt, 1)' calculates the mean along the first dimension (channels).
avg_ref = mean(data_filt, 1);

% Apply the CAR by subtracting the average reference from each channel.
% This re-references the data so that the mean voltage across all channels at any given time point is zero.
data_car = data_filt - avg_ref;

% Plot the PSD of the data after applying the CAR.
figure;
[~,~,~,~,~] = spectopo(data_car, size(data_filt,2), EEG.srate);

%% -------------------------
% Step 4: Epoching around "rt"
% -------------------------

% Extract event latencies and types from the EEG structure.
event_latencies = [EEG.event.latency];
event_types = {EEG.event.type};

% Find the indices of events with the type 'rt' (reaction time).
target_idx = find(strcmp(event_types,'rt'));

% Define the epoch window in samples.
% 'round([-0.2 0.8]*EEG.srate)' calculates the number of samples corresponding to the time window of -200ms to +800ms around the event.
epoch_window = round([-0.2 0.8]*EEG.srate);

% Calculate the length of a single epoch in samples.
epoch_len = diff(epoch_window)+1;

% Preallocate a 3D matrix to store the epochs to improve performance.
% The dimensions will be: channels x epoch length x number of target events.
epochs = nan(size(data_car,1), epoch_len, length(target_idx));

% Loop through each 'rt' event to extract the corresponding epoch.
for i = 1:length(target_idx)
    % Find the center latency of the current event in samples.
    center = round(event_latencies(target_idx(i)));

    % Define the indices of the data points for the current epoch.
    idx = center+epoch_window(1):center+epoch_window(2);

    % Check if the epoch indices are within the bounds of the data matrix to prevent errors.
    if idx(1)>0 && idx(end)<=size(data_car,2)
        % Extract the epoch data and store it in the preallocated matrix.
        epochs(:,:,i) = data_car(:,idx);
    end
end

%% -------------------------
% Step 5: Baseline correction
% -------------------------

% Define the baseline window in samples (0 to 200 ms relative to epoch start).
baseline_idx = 1:round(0.2*EEG.srate);

% Calculate the mean of the baseline period for each channel in each trial.
% 'mean(...,2)' calculates the mean along the second dimension (samples).
% This results in a vector of means for each channel for each trial.
baseline = mean(epochs(:,baseline_idx,:),2);

% Subtract the calculated baseline from each epoch.
% This normalizes the data so that the mean of the baseline period is zero.
% MATLAB's broadcasting handles the subtraction correctly.
epochs_bc = epochs - baseline;

% Extract the baseline data separately for an alternative correction method.
baseline_data = epochs(:,baseline_idx,:);

% Calculate a global baseline by averaging across all channels and all trials within the baseline period.
global_baseline = mean(baseline_data(:));

% Subtract the global baseline from each epoch.
% This is an alternative to the per-channel baseline correction.
epochs_bc_global = epochs - global_baseline;

% -------------------------------------------------------------
% The following section seems to be a separate analysis (spectrogram)
% and not a direct part of the preprocessing pipeline for artifact rejection.
% -------------------------------------------------------------

% Define parameters for the spectrogram analysis.
fs    = EEG.srate;                      % Sampling rate.
win   = hamming(round(0.5*fs));         % Define a 500 ms Hamming window for the spectrogram.
nover = round(0.4*fs);                  % Define a 400 ms overlap between windows.
nfft  = 2^nextpow2(length(win));        % Calculate the FFT length, which must be a power of 2 for efficiency.
nChans  =32;
nTrials =74;

% Preallocate a 4D matrix to store the spectrogram data for all channels and trials.
P_all = [];

% Loop through each channel and each trial to compute the spectrogram.
for ch = 1:nChans
    for tr = 1:nTrials
        % Extract a single trial's data for the current channel.
        sig = squeeze(epochs_bc_global(ch,:,tr));
        
        % Compute the spectrogram.
        % 'spectrogram' returns the power spectral density (P) over time (T) and frequency (F).
        [~,F,T,P] = spectrogram(sig, win, nover, nfft, fs);
        
        % Initialize the preallocated matrix 'P_all' on the first iteration.
        if isempty(P_all)
            P_all = zeros([size(P), nChans, nTrials]);
        end
        
        % Store the computed spectrogram for the current channel and trial.
        P_all(:,:,ch,tr) = P;
    end
end

% === Average across trials and channels ===
% Compute the mean spectrogram by averaging across the channel (3rd) and trial (4th) dimensions.
P_mean = mean(P_all, [3 4]);

% === Plot ===
figure;
% Create a 3D surface plot of the mean spectrogram.
% The x-axis is time, y-axis is frequency, and z-axis is power.
% '10*log10(P_mean)' converts the power to a decibel scale.
surf(T, F, 10*log10(P_mean), 'EdgeColor', 'none');
axis tight; % Adjusts the axes limits to fit the data.
view(0,90);  % Sets the viewing angle to be directly from the top.
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Trial- and Channel-averaged Spectrogram');
colorbar; % Adds a color bar to show the scale of the power values.

%% -------------------------
% Step 6: Improved Artifact Rejection
% -------------------------

% (a) Reject bad trials
% Calculate the variance of each channel within each trial.
% 'var(..., 2)' computes the variance along the second dimension (time).
trial_var = squeeze(var(epochs_bc,0,2));

% Identify trials to reject.
% This finds trials where the mean variance (averaged across channels) is more than 3 standard deviations
% above the global mean variance (across all channels and trials).
trial_reject = (mean(trial_var,1) > mean(mean(trial_var))+3*std(mean(trial_var)));

% Create a logical index for the good trials.
good_trials = ~trial_reject;

% Select only the good trials from the epoch data.
epochs_clean = epochs_bc(:,:,good_trials);

% Print a summary of how many trials were rejected.
fprintf('Rejected %d/%d trials\n', sum(trial_reject), length(trial_reject));

% (b) Reject bad channels (flat or high variance across trials)
% Calculate the variance for each channel across all epochs.
% 'var(..., [2 3])' computes the variance across the time and trial dimensions.
chan_var = squeeze(var(epochs_clean,0,[2 3]));

% Identify channels to reject.
% This finds channels that are either 'flat' (variance less than a tiny value, 1e-6) or
% have a variance that is more than 3 standard deviations above the mean channel variance.
bad_chans = (chan_var < 1e-6) | (chan_var > mean(chan_var)+3*std(chan_var));

% Find the indices of the good channels.
good_chans = find(~bad_chans);

% Select only the good channels from the cleaned epoch data.
epochs_clean = epochs_clean(good_chans,:,:);

% Print a summary of how many channels were rejected.
fprintf('Rejected %d/%d channels\n', sum(bad_chans), length(chan_labels));

% Update the channel labels to reflect the channels that were kept.
chan_labels_clean = chan_labels(good_chans);

%%eegplot
% This section is for plotting and visualization, not part of the core preprocessing.
% It's using the raw EEG data, not the cleaned data.
% Create a vector of channel indices to plot.
channel_indices_to_plot = 1:32;

% Plot the selected channels of the raw, continuous EEG data using 'eegplot'.
% 'srate' specifies the sampling rate.
% 'winlength' sets the duration of the visible window in seconds.
% 'title' provides a descriptive title for the plot.
eegplot(EEG.data(channel_indices_to_plot, :), ...
        'srate', EEG.srate, ...
        'winlength', 10, ...
        'title', ['EEG Data for Task: ' 'example plot' ' (Channels 1-32)']);
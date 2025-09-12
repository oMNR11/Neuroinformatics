%% Loading the BIDS file 

bidsRoot = 'C:\Users\rahma\Downloads\hbn_bids_R1';
bids_demo = bids.layout(bidsRoot,'tolerant',true,'verbose',false);

%% Extracting each individual subject

allSubjectName = {bids_demo.subjects.name};

%% plotting for a single task


% Prepare to plot for the first subject
subject_to_explore = allSubjectName{1};
task_to_explore = 'surroundSupp';

f = bids.query(bids_demo, 'data', ...
                'sub', subject_to_explore, ...
                'task',task_to_explore , ...
                'modality', 'eeg', ...
                'suffix', 'eeg');
%% Load and visualize the EEG data with EEGLAB

% Get the full path to the first file found in the query.
eegFilePath = f{1};


% This will load the data into an EEGLAB structure named 'EEG'.
EEG = pop_loadset(eegFilePath);
data = EEG.data;
chan_labels = {EEG.chanlocs.labels};
figure;
[~,~,~,~,~] = spectopo(data, size(data,2), EEG.srate);


% Plot a segment of the raw EEG data.

channel_indices_to_plot = 1:32;

% Plot only the selected channels with a tighter time window and custom scale
eegplot(EEG.data(channel_indices_to_plot, :), ...
        'srate', EEG.srate, ...
        'winlength', 10, ...
        'title', ['EEG Data for Task: ' task_to_explore ' (Channels 1-32)']);



% %% Loop through each task and plot the EEG data
% 
% 
% % Selecting the first subject
% subject_to_explore = allSubjectName{1};
% 
% % Get all the available tasks for this subject from the BIDS structure
% available_tasks = bids.query(bids_demo, 'tasks', 'sub', subject_to_explore, 'modality', 'eeg');
% 
% if isempty(available_tasks)
%     error(['No EEG tasks found for subject ' subject_to_explore '.']);
% end
% 
% % Loop through each task and generate a plot
% for i = 1:length(available_tasks)
%     current_task = available_tasks{i};
% 
% 
% 
%     % Query for the EEG file for the current task
%     f = bids.query(bids_demo, 'data', ...
%                     'sub', subject_to_explore, ...
%                     'task', current_task, ...
%                     'modality', 'eeg', ...
%                     'suffix', 'eeg');
% 
%     if isempty(f)
%         disp(['No .eeg.set file found for task ' current_task '. Skipping...']);
%         continue;
%     end
% 
%     % Get the full path to the first file found in the query.
%     eegFilePath = f{1};
% 
%     % pop_readbids from EEGLAB loads the BIDS-formatted EEG file.
%     EEG = pop_loadset(eegFilePath);
% 
%     % Plot a segment of the raw EEG data.
%     % We use a new figure for each plot to keep them separate.
%     figure('Name', ['EEG Plot for Task: ' current_task]);
%     eegplot(EEG.data, 'srate', EEG.srate, 'winlength', 5, 'title', ['EEG Data for Task: ' current_task]);
% 
% end


%%
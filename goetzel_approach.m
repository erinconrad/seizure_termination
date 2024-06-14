%{
- thought process:
  - find stim onset and offset using stim channels
  - find ADs using surrounding channels

pseudocode:
- Assume we know what channels we're stimming on (I hope we can get this
from the annotations?)
- Request data every t_stride ms
- For each data pull, look at t_window ms of data
- Look at relative 50 Hz power in the stim channel
- Once it is above threshold, say stim has begun, and start looking for
drop
- Once it falls below new threshold, say stim has stopped, and start
looking at surrounding channels for ADs
- If LL on any surrounding channels exceeds threshold, say AD done
- Reset once we get the next annotation for stim
%}

%% Parameters
f_target = 50;



% rel power in 50 Hz must be above this in both stim channels to say stim
% has begun
stim_start_thresh = 0.1; 

% rel power in 50 Hz must drop below this in both stim channels to say stim
% has stopped
stim_end_thresh = 0.02;

% LL must pass above this in either surround channel to say ADs have
% started
ad_start_thresh = 2e3; 

% LL must go below this in corresponding surround channel to say Ads have
% ended
ad_end_thresh = 1.5e3;

% Consider measuring pre-stim baseline LL and doing relative threshold!!!

%% Main
file_name = 'HUP218_phaseII_D02';
start_time = 193591.76;
end_time = start_time+200;
stim_pair = 'RI04-RI05';
% example of afterdischarges with stim at RI4-5

%% File locs and set path
locations = seizure_termination_paths;
% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
% ieeg stuff
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

%% Get ieeg data
data = download_ieeg_data(file_name,login_name,pwfile,[start_time end_time],1); % 1 means get lots of data
values = data.values;
chLabels = data.chLabels(:,1);
fs = data.fs;
data.start_time = start_time;
data.end_time = end_time;
data.file_name = file_name;


time_window = 0.3; % 100 ms time windows
time_stride = 0.1; % 100 ms stride length

% Define data times
data_stream = values;
data_times = [1:size(values,1)]/fs;
window_length = round(time_window * fs); 
stride = round(time_stride * fs);


%% Get stim channels and surrounding channels
[chs_of_interest,stim_idx] = return_chs_of_interest(stim_pair);
nchs = length(chs_of_interest);
stim_chs = chs_of_interest(logical(stim_idx));
surround_chs = chs_of_interest(~logical(stim_idx));

%% Initialize Goertzel parameters
k = round(f_target * window_length / fs);
omega = 2 * pi * k / window_length;
coeff = 2 * cos(omega);
s_prev = 0;
s_prev2 = 0;
steps = 1:window_length:size(data_stream,1) - window_length + 1;
nsteps = length(steps);
relative_powers = nan(nsteps,nchs);
ll = nan(nsteps,nchs);
stim_on = zeros(nsteps,nchs);
stim_off = zeros(nsteps,nchs);
times = steps/fs;

%% Goertzel algorithm
% Process the signal in windows
count = 0;
for start_idx = 1:window_length:size(data_stream,1) - window_length + 1
    % Extract the window of data
    window_data = data_stream(start_idx:start_idx + window_length - 1,:);
    count = count + 1;

    % notch and bandpass
    [b, a] = butter(2, [58 62]/(fs/2), 'stop');
    window_data = filter(b, a, window_data);
    
    [b, a] = butter(2, [1 70]/(fs/2), 'bandpass');
    window_data = filter(b, a, window_data);


    for j = 1:nchs
        ch = strcmp(chLabels,chs_of_interest{j});
        if sum(ch) == 0,continue; end

        % Apply Goertzel algorithm
        s_prev = 0;
        s_prev2 = 0;
        for n = 1:window_length
            s = window_data(n,ch) + coeff * s_prev - s_prev2;
            s_prev2 = s_prev;
            s_prev = s;
        end
    
        % Calculate power at 50 Hz
        power_50Hz = s_prev2^2 + s_prev^2 - coeff * s_prev * s_prev2;
    
        % Calculate total power in the window
        total_power = sum(window_data(:,ch) .^ 2);
    
        % Calculate relative power
        relative_power = power_50Hz / total_power;
        relative_powers(count,j) = relative_power;

        % ll
        ll(count,j) = sum(abs(diff(window_data(:,ch))));

    end
end

figure
t=tiledlayout(4,3,'tilespacing','tight','padding','tight');
for j =1:nchs
    ch = strcmp(chLabels,chs_of_interest{j});
    nexttile
    plot(data_times,data_stream(:,ch))
    xlim([0 max(data_times)])
    %xlabel('time (s)')
    xticklabels([])
    ylabel('Raw signal')

    nexttile
    plot(steps/fs,relative_powers(:,j))
    hold on
   
    xlim([0 max(data_times)])
    %xlabel('time (s)')
    xticklabels([])
    ylabel('50 Hz relative power')

    nexttile
    plot(steps/fs,ll(:,j))
    xlim([0 max(data_times)])
    %xlabel('time (s)')
    xticklabels([])
    ylabel('LL')
end

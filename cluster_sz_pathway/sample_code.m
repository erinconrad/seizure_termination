%% Notes
%{
times tested: 
[8054.90 8189.90] threshold 2e3 avg_window_sec 2
[10126.72 10201.72] same
[11772.85 11907.85] same
[20634.88 20754.88] ok, a little over sensitive???
[26129.18 26204.18] good
%}

%% Parameters
window_duration = 0.1;         % 100 ms
threshold = 2e3;     
avg_window_sec = 2;
filename = 'HUP247_phaseII';
run_times = [26129.18 26204.18];
sz_ch = {'RP02','RP03'};



%% Paths
locations = seizure_termination_paths;
addpath(genpath(locations.annotation_grabber_folder))
data_folder = [locations.main_folder,'data/'];
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'detections/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% ieeg stuff
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;


%% Load ieeg data
data = download_ieeg_data_sz(filename, login_name, pwfile, run_times, 1);
duration = run_times(2)-run_times(1);
fs = data.fs;
values = data.values;
aT = data.aT;
chLabels = data.chLabels(:,1);
sz_values = values(:,strcmp(chLabels,sz_ch{1}))-values(:,strcmp(chLabels,sz_ch{2})); % bipolar ref
window_size = round(fs * window_duration);
avg_window_samples = avg_window_sec / 0.1;
nsamples = size(values,1);
nwindows = ceil(nsamples/window_size);
times = linspace(0,duration,nsamples);
window_times = linspace(0,duration,nwindows);


%% Test
line_length_history = []; % to store line length of 100 ms windows
seizure_detected = false;
seizure_idx = nan;
ll = nan(nwindows,1);
count = 0;

for start_idx = 1:window_size:(nsamples - window_size + 1)
    count = count + 1;
    segment = sz_values(start_idx:start_idx+window_size-1);

    % Calculate line length for each channel
    line_length = sum(abs(diff(segment)));  

    ll(count) = line_length;
    line_length_history = [line_length_history, line_length];

    % Keep only the last avg_window_samples entries (2 seconds)
    if length(line_length_history) > avg_window_samples
        line_length_history = line_length_history(end - avg_window_samples + 1:end);
    end

    % Compute moving average of line length
    if ~seizure_detected
        if length(line_length_history) == avg_window_samples
            moving_avg = mean(line_length_history);
    
            % Check threshold
            if moving_avg > threshold
                seizure_detected = true;
                seizure_idx = start_idx;
            end
        end
    end

    
end

sz_time = seizure_idx/fs;

figure
tiledlayout(2,1)
nexttile
plot(times,sz_values)
hold on
plot([sz_time sz_time],ylim)

nexttile
plot(window_times,ll)
hold on
plot([sz_time sz_time],ylim)




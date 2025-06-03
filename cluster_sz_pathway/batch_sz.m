%% Parameters
window_duration = 0.1;  % 100 ms
threshold = 2e3;        % Adjust as needed
avg_window_sec = 2;     % Moving average window
filename = 'HUP247_phaseII';
start_time = 59602.5-100;  % Start of your analysis interval
end_time = 76544.3+100;    % End of your analysis interval (example)

chunk_duration = 5 * 60; % 5 minutes in seconds
cooldown_sec = 120;       % 60-second cooldown after a detection
sz_ch = {'RP02','RP03'};

%% Paths and environment setup
locations = seizure_termination_paths;
addpath(genpath(locations.annotation_grabber_folder))
data_folder = [locations.main_folder, 'data/'];
results_folder = [locations.main_folder, 'results/'];
out_folder = [results_folder, 'cluster/'];

scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

%% Initialize detection log
detection_times = [];

%% Outer loop over 5-minute chunks
current_time = start_time;
while current_time < end_time
    % Determine end time for this chunk
    chunk_end_time = min(current_time + chunk_duration, end_time);
    fprintf('Processing chunk: %.2f to %.2f seconds\n', current_time, chunk_end_time);
    
    % Download data chunk
    data = download_ieeg_data_sz(filename, login_name, pwfile, [current_time, chunk_end_time], 1);
    fs = data.fs;
    values = data.values;
    chLabels = data.chLabels(:,1);
    sz_values = values(:, strcmp(chLabels, sz_ch{1})) - values(:, strcmp(chLabels, sz_ch{2}));
    duration = chunk_end_time - current_time;
    
    % Initialize variables for this chunk
    window_size = round(fs * window_duration);
    avg_window_samples = avg_window_sec / window_duration;
    nsamples = size(values, 1);
    nwindows = ceil(nsamples / window_size);
    times = linspace(0, duration, nsamples);
    window_times = linspace(0, duration, nwindows);

    line_length_history = [];
    seizure_detected = false;
    ll = nan(nwindows, 1);
    count = 0;
    
    % Inner loop: analyze 100 ms windows (with cooldown)
    start_idx = 1;
    while start_idx <= (nsamples - window_size + 1)
        count = count + 1;
        segment = sz_values(start_idx:start_idx+window_size-1);
        line_length = sum(abs(diff(segment)));
        ll(count) = line_length;
        line_length_history = [line_length_history, line_length];

        % Keep rolling window of last 2 seconds
        if length(line_length_history) > avg_window_samples
            line_length_history = line_length_history(end - avg_window_samples + 1:end);
        end

        % Check detection
        if ~seizure_detected && length(line_length_history) == avg_window_samples
            moving_avg = mean(line_length_history);
            if moving_avg > threshold
                seizure_detected = true;
                seizure_idx = start_idx;
                seizure_time = current_time + (seizure_idx / fs);

                % Log detection
                detection_times = [detection_times; seizure_time];
                fprintf('Seizure detected at %.2f seconds (absolute)\n', seizure_time);

                % Apply cooldown by skipping samples
                start_idx = start_idx + cooldown_sec * fs;
                % Reset history
                line_length_history = [];
                seizure_detected = false;
                continue; % skip increment at end of loop
            end
        end

        % Increment by window size
        start_idx = start_idx + window_size;
    end
    
    % Move to next chunk
    current_time = chunk_end_time;
end

%% Save detection times to CSV
output_file = fullfile(out_folder, [filename '_detections.csv']);
detection_table = table(detection_times, 'VariableNames', {'SeizureTime_sec'});
writetable(detection_table, output_file);

fprintf('Detections saved to: %s\n', output_file);

%% Optional: plot final detections
if 0
figure
plot(detection_times, 1:length(detection_times), 'o-')
xlabel('Time (s)')
ylabel('Detection number')
title('Seizure detection times')
grid on
end

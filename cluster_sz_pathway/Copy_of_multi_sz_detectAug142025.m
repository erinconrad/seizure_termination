%% measure_accuracy_multi_from_xlsx.m
% Run line-length detector using intervals and channels from an Excel spreadsheet

%% 0. Load spreadsheet -----------------------------------------------------
xls_file = '../../data/cluster_sz_data.xlsx';  
T = readtable(xls_file);

% Ensure filename column is categorical -> string
T.filename = string(T.filename);
T.bipolar_ch1 = string(T.bipolar_ch1);
T.bipolar_ch2 = string(T.bipolar_ch2);

% Unique filenames to loop through
all_filenames = unique(T.filename);

window_duration = 1;       % -> 
avg_window_sec  = 4;         % moving-average window -
chunk_duration  = 5*60;      % 5 min chunks
cooldown_sec    = 180;       % cooldown period in seconds
rel_threshold   = 10;         % # SDs above mean to call sz -> Erin changed from 13 to 10 8/9
notchQ          = 10;        % notch filter
f0              = 60;        % notch filter

%% 1. Paths / env (unchanged) ---------------------------------------------
locations = seizure_termination_paths;
addpath(genpath(locations.annotation_grabber_folder))
data_folder   = [locations.main_folder, 'data/'];
results_folder= [locations.main_folder, 'results/'];
out_folder    = [results_folder, 'cluster/'];

addpath(genpath(locations.script_folder));
addpath(genpath(locations.ieeg_folder));
pwfile     = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

%% 2. Loop over each file --------------------------------------------------
for f = 1:numel(all_filenames)
    filename = all_filenames(f);
    out_csv = fullfile(out_folder, sprintf('%s_detections.csv', filename));

    % Skip this file if already completed
    if isfile(out_csv)
        fprintf('\n✅ Skipping %s (already completed)\n', filename);
        continue
    end
    rows = T(strcmp(T.filename, filename), :);  % all rows for this file

    intervals = [rows.start_time, rows.end_time];
    ch1 = rows.bipolar_ch1;
    ch2 = rows.bipolar_ch2;

    % Estimate threshold using first interval
    baseline_start = intervals(1,1);
    baseline_end   = baseline_start + 300;
    fprintf('\n=== File %s ===\n', filename);
    fprintf('Estimating baseline from %.2f–%.2f s...\n', baseline_start, baseline_end);

    data = download_ieeg_data_sz(filename, login_name, pwfile, [baseline_start, baseline_end], 1);
    fs        = data.fs;
    values    = data.values;
    chLabels  = data.chLabels(:,1);
    sz_values = values(:,strcmp(chLabels, ch1(1))) - values(:,strcmp(chLabels, ch2(1)));


    % Design notch filter
    wo = f0 / (fs/2);     % Normalize frequency
    bw = wo / notchQ;          % Bandwidth
    [b, a] = iirnotch(wo, bw);  % Design notch filter

    % Apply notch to baseline
    sz_values = filtfilt(b, a, sz_values);

    window_size = round(fs * window_duration);
    n_windows   = floor(size(sz_values,1) / window_size);
    ll_vals     = nan(n_windows, 1);

    for w = 1:n_windows
        idx_start = (w-1)*window_size + 1;
        idx_end   = idx_start + window_size - 1;
        if idx_end > size(sz_values,1), break; end
        segment = sz_values(idx_start:idx_end);
        ll_vals(w) = sum(abs(diff(segment)));
    end

    mu_ll  = mean(ll_vals, 'omitnan');
    sigma_ll = std(ll_vals, 'omitnan');
    threshold = mu_ll + rel_threshold * sigma_ll;
    fprintf('Computed threshold = %.1f (mean %.1f + %.1f×std %.1f)\n', ...
        threshold, mu_ll, rel_threshold, sigma_ll);

    %% Detect seizures in each interval
    detection_times_all = [];
    for k = 1:height(rows)
        start_time = rows.start_time(k);
        end_time   = rows.end_time(k);
        szPair     = {ch1(k), ch2(k)};
        fprintf(' Interval %d of %d: %.2f – %.2f s (Channels: %s - %s)\n', ...
                k, height(rows), start_time, end_time, szPair{1}, szPair{2});

        det_times = run_detector_on_interval(filename, start_time, end_time, ...
                       threshold, window_duration, avg_window_sec, ...
                       chunk_duration, cooldown_sec, szPair, login_name, pwfile,...
                       b,a);
        detection_times_all = [detection_times_all; det_times]; %#ok<AGROW>
    end

    %% Save detections
    out_csv = fullfile(out_folder, sprintf('%s_detections.csv', filename));
    writetable(table(detection_times_all, 'VariableNames',{'SeizureTime_sec'}), ...
               out_csv);
    fprintf('  ➜ Saved %d detections to %s\n', numel(detection_times_all), out_csv);
end

function detection_times = run_detector_on_interval(filename, start_time, end_time, ...
                         threshold, window_duration, avg_window_sec, ...
                         chunk_duration, cooldown_sec, sz_ch, login_name, pwfile, ...
                         b,a)

    detection_times = [];
    current_time = start_time;
    last_detection_time = -Inf;
    chunk_counter = 0;

    while current_time < end_time
        chunk_end_time = min(current_time + chunk_duration, end_time);

        % Request garbage collection every 100 chunks
        chunk_counter = chunk_counter + 1;
        if mod(chunk_counter, 100) == 0
            fprintf('🧹 Forcing Java garbage collection...\n');
            java.lang.System.gc();
        end

        data  = download_ieeg_data_sz(filename, login_name, pwfile, ...
                                      [current_time, chunk_end_time], 1);
        fs        = data.fs;
        values    = data.values;
        chLabels  = data.chLabels(:,1);
        sz_values = values(:,strcmp(chLabels,sz_ch{1})) ...
                  - values(:,strcmp(chLabels,sz_ch{2}));

        % Handle NaNs and apply notch
        sz_values(isnan(sz_values)) = nanmean(sz_values);
        if ~any(isnan(sz_values))
            sz_values = filtfilt(b, a, sz_values);
        end

        window_size         = round(fs * window_duration);         % samples in 1s
        n_needed            = round(avg_window_sec / window_duration); % should be 4
        nsamples            = size(values,1);
        line_hist           = [];  % store last n_needed line lengths
        start_idx           = 1;

        while start_idx <= (nsamples - window_size + 1)
            % Compute line length for current 1s window
            segment   = sz_values(start_idx:start_idx+window_size-1);
            line_len  = sum(abs(diff(segment)));
            line_hist = [line_hist, line_len]; %#ok<AGROW>

            % Keep only the last n_needed windows
            if numel(line_hist) > n_needed
                line_hist = line_hist(end-n_needed+1:end);
            end

            % Absolute time aligned to the start of the *current* 1s window
            abs_time_current_win = current_time + (start_idx / fs);

            % Detection rule: all 1s windows in the last 4s must exceed threshold
            if numel(line_hist) == n_needed
                if all(line_hist > threshold)
                    % Timestamp at the start of the 4-window block
                    first_block_start_idx = start_idx - (n_needed-1)*window_size;
                    detection_time = current_time + (first_block_start_idx / fs);

                    if (detection_time - last_detection_time) >= cooldown_sec
                        detection_times(end+1,1) = detection_time; %#ok<AGROW>
                        fprintf('    Detected at %.2f s (start of 4-window block)\n', detection_time);
                        last_detection_time = detection_time;

                        % Jump ahead by cooldown to avoid re-detecting the same event
                        start_idx = start_idx + round(cooldown_sec * fs);
                        line_hist = [];
                        continue
                    end
                end
            end

            % Advance by one non-overlapping 1s window
            start_idx = start_idx + window_size;
        end

        current_time = chunk_end_time;
    end
end

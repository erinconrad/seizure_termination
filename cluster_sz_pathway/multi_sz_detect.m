%% measure_accuracy_multi.m
% Run line‑length detector over many files/intervals, each with its own threshold

%% 0. Config: list the files, intervals, and thresholds --------------------
% Each row: filename, cell array of [start end] rows, threshold
fileSpecs = struct( ...
    'filename',  {'HUP247_phaseII',...
                'HUP251_phaseII',...
                'HUP266_phaseII'}, ...
    'intervals', { [60052.43-400, 80786.43+100],...
                [30043.37207-400, 35440.62988+100],...
                [253184.72-400, 301581.35+100], ...
                   }, ...
    'sz_ch',     {{'RP02','RP03'},...
                {'LB01','LB02'}
                {'RC01','RC02'}} ...
);

window_duration = 0.1;       % 100 ms
avg_window_sec  = 2;         % moving‑average window
chunk_duration  = 5*60;      % 5 min chunks
cooldown_sec    = 120;       % cooldown

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

%% 2. Loop over files ------------------------------------------------------
for f = 1:numel(fileSpecs)
    spec        = fileSpecs(f);
    filename    = spec.filename;
    intervals   = spec.intervals;
    szPair      = spec.sz_ch;

    %% NEW: Estimate threshold from first 5 minutes
    baseline_start = intervals(1,1);         % assume this is safe pre-seizure
    baseline_end   = baseline_start + 300;   % 5 min = 300 sec

    fprintf('\n=== File %s ===\n', filename);
    fprintf('Estimating baseline from %.2f–%.2f s...\n', baseline_start, baseline_end);

    % Get baseline data
    data = download_ieeg_data_sz(filename, login_name, pwfile, [baseline_start, baseline_end], 1);
    fs        = data.fs;
    values    = data.values;
    chLabels  = data.chLabels(:,1);
    sz_values = values(:,strcmp(chLabels,szPair{1})) - values(:,strcmp(chLabels,szPair{2}));

    % Compute line length per window
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
    threshold = mu_ll + 4 * sigma_ll;    % Tune the multiplier if needed

    fprintf('Computed threshold = %.1f (mean %.1f + 5×std %.1f)\n', threshold, mu_ll, sigma_ll);

    %% Detect seizures in intervals
    detection_times_all = [];

    for k = 1:size(intervals,1)
        start_time = intervals(k,1);
        end_time   = intervals(k,2);
        fprintf(' Interval %d of %d: %.2f – %.2f s\n', ...
                k, size(intervals,1), start_time, end_time);
        
        det_times = run_detector_on_interval(filename, start_time, end_time, ...
                       threshold, window_duration, avg_window_sec, ...
                       chunk_duration, cooldown_sec, szPair, login_name, pwfile);
        detection_times_all = [detection_times_all; det_times]; %#ok<AGROW>
    end

    %% Save detections
    out_csv = fullfile(out_folder, sprintf('%s_detections.csv', filename));
    writetable(table(detection_times_all, 'VariableNames',{'SeizureTime_sec'}), ...
               out_csv);
    fprintf('  ➜ Saved %d detections to %s\n', numel(detection_times_all), out_csv);
end


%% ------------------------------------------------------------------------
% Helper: run detector over one [start end] interval
function detection_times = run_detector_on_interval(filename, start_time, end_time, ...
                         threshold, window_duration, avg_window_sec, ...
                         chunk_duration, cooldown_sec, sz_ch, login_name, pwfile)

    detection_times = [];
    current_time = start_time;

    while current_time < end_time
        % --- fetch data chunk -------------------------------------------
        chunk_end_time = min(current_time + chunk_duration, end_time);
        data  = download_ieeg_data_sz(filename, login_name, pwfile, ...
                                      [current_time, chunk_end_time], 1);
        fs        = data.fs;
        values    = data.values;
        chLabels  = data.chLabels(:,1);
        sz_values = values(:,strcmp(chLabels,sz_ch{1})) ...
                  - values(:,strcmp(chLabels,sz_ch{2}));

        % --- window prep -------------------------------------------------
        window_size        = round(fs * window_duration);
        avg_window_samples = avg_window_sec / window_duration;
        nsamples           = size(values,1);
        seizure_detected   = false;
        line_hist          = [];
        start_idx          = 1;

        %% --- slide 100 ms windows --------------------------------------
        while start_idx <= (nsamples - window_size + 1)
            segment   = sz_values(start_idx:start_idx+window_size-1);
            line_len  = sum(abs(diff(segment)));
            line_hist = [line_hist, line_len];      %#ok<AGROW>

            if numel(line_hist) > avg_window_samples
                line_hist = line_hist(end-avg_window_samples+1:end);
            end

            if ~seizure_detected && numel(line_hist) == avg_window_samples
                if mean(line_hist) > threshold
                    seizure_detected = true;
                    seizure_time = current_time + (start_idx / fs);
                    detection_times(end+1,1) = seizure_time; %#ok<AGROW>
                    fprintf('    Detected at %.2f s (abs)\n', seizure_time);

                    start_idx = start_idx + cooldown_sec * fs;  % cooldown
                    line_hist = [];  seizure_detected = false;
                    continue
                end
            end
            start_idx = start_idx + window_size;
        end
        current_time = chunk_end_time;  % next chunk
    end
end

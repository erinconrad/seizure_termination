%% measure_accuracy_multi_from_xlsx_adaptive.m
% Adaptive baseline + robust z-score + 4×1s persistence
% Notes:
%   - Baseline tracks the *feature* (line length per 1s window), not the raw signal
%   - Dual EWMA trackers (fast/slow) for center and scale (mean abs deviation proxy)
%   - Updates are gated to avoid contaminating baseline with candidate/artifact windows
%   - Detection: ALL of the last 5 one-second windows must exceed T_high (in z-space)
%   - Hysteresis: event ends when z falls below T_low

%% 0. Load spreadsheet -----------------------------------------------------
xls_file = '../../data/cluster_sz_data.xlsx';
T = readtable(xls_file);

% Ensure filename column is categorical -> string
T.filename    = string(T.filename);
T.bipolar_ch1 = string(T.bipolar_ch1);
T.bipolar_ch2 = string(T.bipolar_ch2);

% Unique filenames to loop through
all_filenames = unique(T.filename);

% --- Params you had + a few new ones for the adaptive baseline -----------
window_duration = 1;          % seconds per window
chunk_duration  = 5*60;       % 5 min chunks
cooldown_sec    = 180;        % do not re-detect within 180 s
f0 = 60; NOTCH_Q = 35; harmonics = [1 2];   % 60 & 120 Hz notch
LPF_CUTOFF = 50;                             

% Adaptive baseline params (tweak here)
T_high         = 8;          % enter threshold in z units
T_low          = 6;           % exit threshold in z units (hysteresis)
K_needed       = 5;           % ALL last 5 windows must be > T_high
refractory_sec = 10;          % freeze updates after an event

% EWMA half-lives (convert to per-window alphas below)
tau_fast_s     = 60;          % ~1 min fast tracker half-life
tau_slow_s     = 30*60;       % ~30 min slow tracker half-life
beta_mix       = 0.2;         % combine slow/fast: mu = (1-b)*slow + b*fast

% Gating: treat very large spikes as candidates (don’t update baseline on them)
gate_z_for_update = T_high - 0.5;        % if z>8, skip baseline updates

%% 1. Paths / env (unchanged) ---------------------------------------------
locations = seizure_termination_paths;
addpath(genpath(locations.annotation_grabber_folder))
data_folder    = [locations.main_folder, 'data/'];
results_folder = [locations.main_folder, 'results/'];
out_folder     = [results_folder, 'cluster/'];

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

    % --- Precompute filter objects once per file -------------------------
    % Design notch (60, 120 Hz if below Nyquist)
    data_probe = download_ieeg_data_sz(filename, login_name, pwfile, [intervals(1,1), min(intervals(1,1)+5, intervals(1,2))], 1);
    fs = data_probe.fs;

    sos_notch = [];
    g_notch = 1;
    for h = harmonics
        fr = h*f0;
        if fr < fs/2-1
            wo = fr/(fs/2);
            bw = wo/NOTCH_Q;
            [b1,a1] = iirnotch(wo,bw);
            [sos, g] = tf2sos(b1,a1);
            sos_notch = [sos_notch; sos]; %#ok<AGROW>
            g_notch = g_notch * g;
        end
    end
    % 4th-order LPF
    [blp,alp] = butter(4, LPF_CUTOFF/(fs/2), 'low');
    [sos_lpf, g_lpf] = tf2sos(blp, alp);

    % Convert half-lives to per-window alphas
    alpha_fast = 1 - 2^(-window_duration / tau_fast_s);
    alpha_slow = 1 - 2^(-window_duration / tau_slow_s);

    fprintf('\n=== File %s ===\n', filename);

    %% Detect seizures in each interval
    detection_times_all = [];
    for k = 1:height(rows)
        start_time = rows.start_time(k);
        end_time   = rows.end_time(k);
        szPair     = {ch1(k), ch2(k)};
        fprintf(' Interval %d of %d: %.2f – %.2f s (Channels: %s - %s)\n', ...
                k, height(rows), start_time, end_time, szPair{1}, szPair{2});

        % Initialize per-interval adaptive baseline state
        st = init_adaptive_state(alpha_fast, alpha_slow, beta_mix, ...
                                 T_high, T_low, K_needed, refractory_sec, ...
                                 window_duration);

        % Run detector on this interval (returns updated state too)
        [det_times, ~] = run_detector_on_interval_adaptive( ...
                       filename, start_time, end_time, ...
                       window_duration, ...
                       chunk_duration, cooldown_sec, szPair, login_name, pwfile, ...
                       sos_notch, g_notch, sos_lpf, g_lpf, ...
                       st, gate_z_for_update);

        detection_times_all = [detection_times_all; det_times]; %#ok<AGROW>
    end

    %% Save detections
    out_csv = fullfile(out_folder, sprintf('%s_detections.csv', filename));
    writetable(table(detection_times_all, 'VariableNames',{'SeizureTime_sec'}), out_csv);
    fprintf('  ➜ Saved %d detections to %s\n', numel(detection_times_all), out_csv);
end


%% --- Helper: initialize adaptive baseline state -------------------------
function st = init_adaptive_state(alpha_fast, alpha_slow, beta_mix, T_high, T_low, K_needed, refractory_sec, win_sec)
    st.mu_fast = NaN;               % will be bootstrapped on first window
    st.d_fast  = NaN;
    st.mu_slow = NaN;
    st.d_slow  = NaN;

    st.alpha_fast = alpha_fast;
    st.alpha_slow = alpha_slow;
    st.beta       = beta_mix;

    st.T_high = T_high;
    st.T_low  = T_low;
    st.K      = K_needed;

    st.refractory_sec = refractory_sec;

    st.in_event = false;
    st.refrac_until_abs_t = -Inf;

    st.buf = false(1, st.K);     % ring buffer of last K above-threshold flags
    st.buf_idx = 1;

    st.last_detection_time = -Inf;
    st.win_sec = win_sec;
end


%% --- Main detector with adaptive baseline -------------------------------
function [detection_times, st] = run_detector_on_interval_adaptive( ...
                         filename, start_time, end_time, ...
                         window_duration,  ...
                         chunk_duration, cooldown_sec, sz_ch, login_name, pwfile, ...
                         sos_notch, g_notch, sos_lpf, g_lpf, ...
                         st, gate_z_for_update)

    detection_times = [];
    current_time = start_time;
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

        % Handle NaNs
        if any(isnan(sz_values))
            nz = ~isnan(sz_values);
            if any(nz)
                sz_values(~nz) = mean(sz_values(nz));
            else
                % all NaN—skip this chunk
                current_time = chunk_end_time;
                continue
            end
        end

        % Apply cascade: notch -> LPF (forward-only, low latency)
        if ~isempty(sos_notch), sz_values = sosfilt(sos_notch, g_notch * sz_values); end
        if ~isempty(sos_lpf),   sz_values = sosfilt(sos_lpf,   g_lpf   * sz_values); end

        window_size   = round(fs * window_duration);   % samples / 1 s window
        nsamples      = size(sz_values,1);
        start_idx     = 1;

        while start_idx <= (nsamples - window_size + 1)
            seg      = sz_values(start_idx:start_idx+window_size-1);
            line_len = sum(abs(diff(seg)));

            % --- Bootstrap baseline on very first window if needed
            if isnan(st.mu_fast)
                st.mu_fast = line_len;
                st.mu_slow = line_len;
                st.d_fast  = max(1, 0.5*abs(line_len));   % conservative nonzero
                st.d_slow  = st.d_fast;
            end

            % Combine baseline center & scale
            mu  = (1 - st.beta)*st.mu_slow + st.beta*st.mu_fast;
            dev = (1 - st.beta)*st.d_slow  + st.beta*st.d_fast + eps;

            % Robust z for this window
            z = (line_len - mu) / dev;

            % --- Detection rule: ALL last K windows above T_high -----------
            is_hi = (z > st.T_high);
            st.buf(st.buf_idx) = is_hi;
            st.buf_idx = 1 + mod(st.buf_idx, st.K);

            abs_time_current_win = current_time + (start_idx / fs);                  % window start
            abs_time_end_of_win  = current_time + ((start_idx + window_size) / fs);  % window end

            if ~st.in_event
                if all(st.buf) && (abs_time_end_of_win - st.last_detection_time) >= cooldown_sec ...
                   && abs_time_end_of_win >= st.refrac_until_abs_t
                    % Mark detection time at the *end* of the 5-window block
                    detection_times(end+1,1) = abs_time_end_of_win; %#ok<AGROW>
                    fprintf('    Detected at %.2f s (end of 5-window block)\n', abs_time_end_of_win);
                    st.in_event = true;
                    st.last_detection_time = abs_time_end_of_win;
                    st.refrac_until_abs_t  = abs_time_end_of_win + st.refractory_sec;

                    % Optional: fast-forward to near refractory boundary to save compute
                    % But keep it simple/real-time friendly; just clear the buffer
                    st.buf(:) = false;
                end
            else
                % End event with hysteresis
                if z < st.T_low
                    st.in_event = false;
                    st.buf(:)   = false;
                    % refractory already set on entry
                end
            end

            % --- Baseline update gating (skip if candidate/artifact) -------
            % Skip updates if currently in event, in refractory, or this window is a large outlier
            can_update = (~st.in_event) && (abs_time_current_win >= st.refrac_until_abs_t) && (z <= gate_z_for_update);

            if can_update
                aF = st.alpha_fast; aS = st.alpha_slow;

                % Update centers
                st.mu_fast = st.mu_fast + aF*(line_len - st.mu_fast);
                st.mu_slow = st.mu_slow + aS*(line_len - st.mu_slow);

                % Update robust deviation proxies (mean absolute deviation)
                st.d_fast  = st.d_fast  + aF*(abs(line_len - st.mu_fast) - st.d_fast);
                st.d_slow  = st.d_slow  + aS*(abs(line_len - st.mu_slow) - st.d_slow);

                % Keep scales positive/nonzero
                st.d_fast = max(st.d_fast, 1e-6);
                st.d_slow = max(st.d_slow, 1e-6);
            end

            % Advance by one non-overlapping 1 s window
            start_idx = start_idx + window_size;
        end

        current_time = chunk_end_time;
    end
end

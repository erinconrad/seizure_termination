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

filename = 'HUP184_phaseII'


function detection_times = run_detector_on_interval_test_inner(filename, start_time, end_time, ...
                         threshold, window_duration, avg_window_sec, ...
                         chunk_duration, cooldown_sec, sz_ch, login_name, pwfile,...
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

        sz_values(isnan(sz_values)) = nanmean(sz_values); % bug wiht nans

        if ~any(isnan(sz_values))
            sz_values = filtfilt(b, a, sz_values); % apply notch
        end
        

        window_size        = round(fs * window_duration);
        avg_window_samples = avg_window_sec / window_duration;
        nsamples           = size(values,1);
        line_hist          = [];
        start_idx          = 1;

        while start_idx <= (nsamples - window_size + 1)
            segment   = sz_values(start_idx:start_idx+window_size-1);
            line_len  = sum(abs(diff(segment)));
            line_hist = [line_hist, line_len]; %#ok<AGROW>

            if numel(line_hist) > avg_window_samples
                line_hist = line_hist(end-avg_window_samples+1:end);
            end

            abs_time = current_time + (start_idx / fs);
            if numel(line_hist) == avg_window_samples
                if mean(line_hist) > threshold && (abs_time - last_detection_time >= cooldown_sec)
                    detection_times(end+1,1) = abs_time; %#ok<AGROW>
                    fprintf('    Detected at %.2f s (abs)\n', abs_time);
                    last_detection_time = abs_time;

                    start_idx = start_idx + round(cooldown_sec * fs);
                    line_hist = [];
                    continue
                end
            end
            start_idx = start_idx + window_size;
        end
        current_time = chunk_end_time;
    end
end
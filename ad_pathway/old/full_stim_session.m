%{
- to dos
  - consider doing ll threshold relative to baseline in that channel
  - eventually need to improve AD off detection
  - test for many more patients
  - talk to folks about benefits/costs of different window sizes and stride
  lengths, what makes sense in continuous scenario
%}

clear

%% Parameters
% Goertzel parameters
f_target = 50;
time_window = 0.2; % 100 ms time windows
time_stride = 0.2; % 100 ms stride length

% rel power in 50 Hz must be above this in both stim channels to say stim
% has begun
stim_start_thresh = 20; 

% rel power in 50 Hz must drop below this in both stim channels to say stim
% has stopped
stim_end_thresh = 5;

% LL must pass above this in either surround channel to say ADs have
% started
ad_start_thresh = 2e3; 

% LL must go below this in corresponding surround channel to say Ads have
% ended
ad_end_thresh = 0.5e3;

stop_looking_time = 5; % if ADs haven't started by now, stop looking

% Consider measuring pre-stim baseline LL and doing relative threshold!!!

%% Main
file_name = 'HUP218_HFS';
start_time = 21238.33;%20736.07; 
end_time = 21238.33+10;%20736.07+15;%23762.34;

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
chLabels = decompose_labels(chLabels);
fs = data.fs;
data.start_time = start_time;
data.end_time = end_time;
data.file_name = file_name;
aT = data.aT;

% Define data times
data_stream = values;
data_times = [1:size(values,1)]/fs;
ieeg_times = linspace(start_time,end_time,size(values,1));
window_length = round(time_window * fs); 
stride = round(time_stride * fs);


%% Start looking for stim and ADs
% Initialize Goertzel parameters
k = round(f_target * window_length / fs);
omega = 2 * pi * k / window_length;
coeff = 2 * cos(omega);
s_prev = 0;
s_prev2 = 0;

% Filters
[b_notch,a_notch] = butter(2, [58 62]/(fs/2), 'stop');
z_notch = zeros(max(length(a_notch),length(b_notch))-1,1); % filter state

[b_pass,a_pass] = butter(2, [1 70]/(fs/2), 'bandpass');
z_pass = zeros(max(length(a_pass),length(b_pass))-1,1); % filter state

% Initialize conditions
stim_begun = zeros(2,1);
stim_ended = zeros(2,1);

ad_begun = 0;
ad_ended = 0;
ad_ch = nan;

% Initialize time stamps
time_stamps = cell(0,3);
stim_pair = [];
start_idx = 1;
first_time = 1;

% Loop over data stream
for i = window_length+1:stride:length(data_stream)

    % Define data window
    curr_values = data_stream(i-window_length:i,:);
    curr_times = [start_time + data_times(i-window_length),...
        start_time + data_times(i)];

    % Check for closed or open relay annotations in this window
    new_stim_pair = find_relay_anns(aT,curr_times);
    old_stim_pair = stim_pair;
    if ~isempty(new_stim_pair)
        stim_pair = new_stim_pair;

    end
        
    if ~isempty(old_stim_pair) && ~strcmp(old_stim_pair,new_stim_pair) && ~isempty(new_stim_pair)
        end_idx = i;
        time_stamp_times = cell2mat(time_stamps(:,1));
        curr_time_stamp_indices = time_stamp_times >= start_idx & time_stamp_times <= end_idx;
        plot_stims(chLabels,chs_of_interest,data_times,data_stream,...
            time_stamps(curr_time_stamp_indices,:),fs,[start_idx end_idx])
        old_stim_pair
        new_stim_pair
        first_time = 0;
        pause
        start_idx = i;
    end

    if isempty(stim_pair)
        continue
    end


    if ~isempty(new_stim_pair) || first_time == 1
        % Get stim channels and surrounding channels
        [chs_of_interest,stim_idx] = return_chs_of_interest(stim_pair);
        nchs = length(chs_of_interest);
        stim_chs = chs_of_interest(logical(stim_idx));
        surround_chs = chs_of_interest(~logical(stim_idx));

        % Get baseline characteristics on channels
    end

    % notch and bandpass
    [curr_values,z_notch] = filter(b_notch, a_notch, curr_values,z_notch);
    [curr_values,z_pass] = filter(b_pass, a_pass, curr_values,z_pass);

    % if stim hasn't begun on both channel
    if any(stim_begun==0) 

        % Look for the beginning of stim

        % Loop over stim channels
        for j = 1:2
            if stim_begun(j) ~=0, continue; end % skip if stim already started on that ch
            ch = strcmp(chLabels,stim_chs{j});
            if sum(ch) == 0,error('cannot find stim ch'); end

            % restrict to channel
            data_ch = curr_values(:,ch);

            % Apply Goertzel algorithm
            s_prev = 0;
            s_prev2 = 0;
            for n = 1:window_length
                s = curr_values(n,ch) + coeff * s_prev - s_prev2;
                s_prev2 = s_prev;
                s_prev = s;
            end

            % Calculate power at 50 Hz
            power_50Hz = s_prev2^2 + s_prev^2 - coeff * s_prev * s_prev2;

            % Calculate total power in the window
            total_power = sum(data_ch .^ 2);
            rel_50 = power_50Hz/total_power;

            if rel_50 > stim_start_thresh
                stim_begun(j) = i; % say stim has begun on that channel
                time_stamps = [time_stamps;{i},stim_chs(j),{'Stim start'}];
            end
        end

    % If stim has begun on both channels and hasn't ended on both channels
    elseif all(stim_begun~=0) && any(stim_ended==0)
        % Look for end of stim
        

        % Loop over stim channels
        for j = 1:2
            if stim_ended(j) ~=0, continue; end % skip if stim already ended on that ch
            ch = strcmp(chLabels,stim_chs{j});
            if sum(ch) == 0,error('cannot find stim ch'); end

            % restrict to channel
            data_ch = curr_values(:,ch);

            % Apply Goertzel algorithm
            s_prev = 0;
            s_prev2 = 0;
            for n = 1:window_length
                s = curr_values(n,ch) + coeff * s_prev - s_prev2;
                s_prev2 = s_prev;
                s_prev = s;
            end

            % Calculate power at 50 Hz
            power_50Hz = s_prev2^2 + s_prev^2 - coeff * s_prev * s_prev2;

            % Calculate total power in the window
            total_power = sum(data_ch .^ 2);
            rel_50 = power_50Hz/total_power;

            if rel_50 < stim_end_thresh
                stim_ended(j) = i; % say stim has ended on that channel
                time_stamps = [time_stamps;{i},stim_chs(j),{'Stim end'}];
            end
        end

    % If stim has begun and ended on both channels and ADs not yet begun on
    % either
    elseif all(stim_begun~=0) && all(stim_ended~=0) && ad_begun == 0

        % Loop over surround channels
        if i/fs - max(stim_ended)/fs < stop_looking_time
            for j = 1:2
                ch = strcmp(chLabels,surround_chs{j});
                if sum(ch) == 0,continue; end
                data_ch = curr_values(:,ch);
                ll = sum(abs(diff(data_ch)));
    
                if ll > ad_start_thresh
                    ad_begun = i; % AD started
                    ad_ch = surround_chs{j};
                    time_stamps = [time_stamps;{i},ad_ch,{'ADs start'}];
    
                end
            end
        end

        % Also loop over stim to see if stim starts before ADs
        for j = 1:2
            ch = strcmp(chLabels,stim_chs{j});
            if sum(ch) == 0,error('cannot find stim ch'); end

            % restrict to channel
            data_ch = curr_values(:,ch);

            % Apply Goertzel algorithm
            s_prev = 0;
            s_prev2 = 0;
            for n = 1:window_length
                s = curr_values(n,ch) + coeff * s_prev - s_prev2;
                s_prev2 = s_prev;
                s_prev = s;
            end

            % Calculate power at 50 Hz
            power_50Hz = s_prev2^2 + s_prev^2 - coeff * s_prev * s_prev2;

            % Calculate total power in the window
            total_power = sum(data_ch .^ 2);
            rel_50 = power_50Hz/total_power;

            if rel_50 > stim_start_thresh
                stim_begun(j) = i; % say stim has begun on that channel
 
                stim_ended(j) = 0; % reset stim stopped to 0
                time_stamps = [time_stamps;{i},stim_chs(j),{'Stim start'}];
            end
        end

    % If stim has begun and ended, and ADs have begun but not yet ended    
    elseif all(stim_begun~=0) && all(stim_ended~=0) && ad_begun ~=0 && ad_ended == 0

        % Look at AD channel
        ch = strcmp(chLabels,ad_ch);

        if sum(ch) == 0,error('cannot find ad ch'); end

        data_ch = curr_values(:,ch);
        ll = sum(abs(diff(data_ch)));

        if ll < ad_end_thresh
            ad_ended = i; % ended ADs
            time_stamps = [time_stamps;{i},ad_ch,{'ADs end'}];
            % Plot stuff
            %plot_stims(chLabels,chs_of_interest,data_times,data_stream,...
  %  stim_begun,stim_ended,ad_begun,ad_ended,ad_ch,fs)
            %error('look')

        end

        % Also loop over stim to see if stim starts before ADs end
        for j = 1:2
            ch = strcmp(chLabels,stim_chs{j});
            if sum(ch) == 0,error('cannot find stim ch'); end

            % restrict to channel
            data_ch = curr_values(:,ch);

            % Apply Goertzel algorithm
            s_prev = 0;
            s_prev2 = 0;
            for n = 1:window_length
                s = curr_values(n,ch) + coeff * s_prev - s_prev2;
                s_prev2 = s_prev;
                s_prev = s;
            end

            % Calculate power at 50 Hz
            power_50Hz = s_prev2^2 + s_prev^2 - coeff * s_prev * s_prev2;

            % Calculate total power in the window
            total_power = sum(data_ch .^ 2);
            rel_50 = power_50Hz/total_power;

            if rel_50 > stim_start_thresh
                stim_begun(j) = i; % say stim has begun on that channel
 
                stim_ended(j) = 0; % reset stim stopped to 0
                ad_begun = 0; % also reset ad
                time_stamps = [time_stamps;{i},stim_chs(j),{'Stim start'}];
            end
        end

    elseif all(stim_begun~=0) && all(stim_ended~=0) && ad_begun ~=0 && ad_ended ~=0

        % everything done! Reset
        stim_begun = zeros(2,1);
        stim_ended = zeros(2,1);
        
        ad_begun = 0;
        ad_ended = 0;
        ad_ch = nan;
    else
        error('what')

    end

end


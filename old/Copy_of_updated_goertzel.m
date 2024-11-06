%{
- to dos
  - consider doing ll threshold relative to baseline in that channel
  - eventually need to improve AD off detection
  - test for many more patients
  - talk to folks about benefits/costs of different window sizes and stride
  lengths, what makes sense in continuous scenario
%}

%% Parameters
% Goertzel parameters
f_target = 50;
time_window = 0.3; % 100 ms time windows
time_stride = 0.1; % 100 ms stride length

% rel power in 50 Hz must be above this in both stim channels to say stim
% has begun
stim_start_thresh = 10; 

% rel power in 50 Hz must drop below this in both stim channels to say stim
% has stopped
stim_end_thresh = 5;

% LL must pass above this in either surround channel to say ADs have
% started
ad_start_thresh = 3e3; 

% LL must go below this in corresponding surround channel to say Ads have
% ended
ad_end_thresh = 2e3;

% Consider measuring pre-stim baseline LL and doing relative threshold!!!

%% Main
file_name = 'HUP218_phaseII_D02';
start_time = 193591.76+0; % add 10 to get to 2nd; add 60 to get to 3rd; 80 4th; 130 5th
end_time = 193591.76+200;
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

%% Start looking for stim and ADs
% Initialize Goertzel parameters
k = round(f_target * window_length / fs);
omega = 2 * pi * k / window_length;
coeff = 2 * cos(omega);
s_prev = 0;
s_prev2 = 0;

% Initialize conditions
stim_begun = zeros(2,1);
stim_ended = zeros(2,1);

ad_begun = 0;
ad_ended = 0;
ad_ch = nan;

% Loop over data stream
for i = window_length+1:stride:length(data_stream)

    % Define data window
    data = data_stream(i-window_length:i,:);

    % notch and bandpass
    [b, a] = butter(2, [58 62]/(fs/2), 'stop');
    data = filter(b, a, data);
    
    [b, a] = butter(2, [1 70]/(fs/2), 'bandpass');
    data = filter(b, a, data);

    % if stim hasn't begun on both channel
    if any(stim_begun==0) 

        % Look for the beginning of stim

        % Loop over stim channels
        for j = 1:2
            ch = strcmp(chLabels,stim_chs{j});
            if sum(ch) == 0,error('cannot find stim ch'); end

            % restrict to channel
            data_ch = data(:,ch);

            % Apply Goertzel algorithm
            s_prev = 0;
            s_prev2 = 0;
            for n = 1:window_length
                s = data(n,ch) + coeff * s_prev - s_prev2;
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
            end
        end

    % If stim has begun on both channels and hasn't ended on both channels
    elseif all(stim_begun~=0) && any(stim_ended==0)
        % Look for end of stim
        

        % Loop over stim channels
        for j = 1:2
            ch = strcmp(chLabels,stim_chs{j});
            if sum(ch) == 0,error('cannot find stim ch'); end

            % restrict to channel
            data_ch = data(:,ch);

            % Apply Goertzel algorithm
            s_prev = 0;
            s_prev2 = 0;
            for n = 1:window_length
                s = data(n,ch) + coeff * s_prev - s_prev2;
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
            end
        end

    % If stim has begun and ended on both channels and ADs not yet begun on
    % either
    elseif all(stim_begun~=0) && all(stim_ended~=0) && ad_begun == 0

        % Loop over surround channels
        for j = 1:2
            ch = strcmp(chLabels,surround_chs{j});
            if sum(ch) == 0,continue; end
            data_ch = data(:,ch);
            ll = sum(abs(diff(data_ch)));

            if ll > ad_start_thresh
                ad_begun = i; % AD started
                ad_ch = surround_chs{j};


            end
        end

    % If stim has begun and ended, and ADs have begun but not yet ended    
    elseif all(stim_begun~=0) && all(stim_ended~=0) && ad_begun ~=0 && ad_ended == 0

        % Look at AD channel
        ch = strcmp(chLabels,ad_ch);

        if sum(ch) == 0,error('cannot find ad ch'); end

        data_ch = data(:,ch);
        ll = sum(abs(diff(data_ch)));

        if ll < ad_end_thresh
            ad_ended = i; % ended ADs

            % Plot stuff
            plot_stims(chLabels,chs_of_interest,data_times,data_stream,...
    stim_begun,stim_ended,ad_begun,ad_ended,ad_ch,fs)
            error('look')

        end

    elseif all(stim_begun~=0) && all(stim_ended~=0) && ad_begun ~=0 && ad_ended ~=0

        % everything done! Reset
    else
        error('what')

    end

end
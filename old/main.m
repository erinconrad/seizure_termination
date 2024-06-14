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
time_window = 0.3; % 100 ms time windows
time_stride = 0.1; % 200 ms stride length

% rel power in 50 Hz must be above this in both stim channels to say stim
% has begun
stim_start_thresh = 0.1; 

% rel power in 50 Hz must drop below this in both stim channels to say stim
% has stopped
stim_end_thresh = 0.02;

% LL must pass above this in either surround channel to say ADs have
% started
ad_start_thresh = 4e3; 

% LL must go below this in corresponding surround channel to say Ads have
% ended
ad_end_thresh = 2e3;

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


%% WHAT
if 1
steps = window_length+1:stride:length(data_stream);
nsteps = length(steps);

rel_50 = nan(nsteps,nchs);
ll = nan(nsteps,nchs);

stim_on = zeros(nsteps,nchs);
stim_off = zeros(nsteps,nchs);
ad_on = zeros(nsteps,nchs);



count = 0;

for i = window_length+1:stride:length(data_stream)
    data = data_stream(i-window_length:i,:);
    count = count + 1;     
    % notch and bandpass
    [b, a] = butter(2, [58 62]/(fs/2), 'stop');
    data = filter(b, a, data);
    
    [b, a] = butter(2, [1 70]/(fs/2), 'bandpass');
    data = filter(b, a, data);

    % 50 hz filter
    [b, a] = butter(2, [48 52]/(fs/2), 'bandpass');
    filtered_50 = filter(b, a, data);
    for j = 1:nchs
        
        ch = strcmp(chLabels,chs_of_interest{j});
        if sum(ch) == 0,continue; end

        data_ch = data(:,ch);
        filtered_50_ch = filtered_50(:,ch);
           
    
        rel_50(count,j) = sum(filtered_50_ch.^2)./sum(data_ch.^2);
        ll(count,j) = sum(abs(diff(data_ch)));
    
        if rel_50(count,j) > thresh_50
            % turn it on
            
            stim_on(count,j) = 1;
                
            
        else
            if count > 1
                % if stim was just on, turn it off
                if stim_on(count-1,j) == 1
                    stim_on(count,j) = 0;
                    stim_off(count,j) = 1;
                else
                    % if stim off, look for ads
                    if stim_off(count-1,j) == 1
                        if ll(count,j) > thresh_ll
                            ad_on(count,j) = 1;
                            %error('look')
                        else
                            ad_on(count,j) = 0;
        
                        end
                    end
                end
            end
    
        end
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
    plot(steps/fs,rel_50(:,j))
    hold on
    %arrayfun(@(x) xline(x, 'Color', 'g', 'LineWidth', 2), steps(stim_on(:,j)==1)/fs);
    %plot(steps(stim_on(:,j)==1)/fs,rel_50(stim_on(:,j)==1,j),'go')
    arrayfun(@(x) xline(x, 'Color', 'r', 'LineWidth', 2), steps(stim_off(:,j)==1)/fs);
    %plot([steps(stim_off(:,j)==1)/fs steps(stim_off(:,j)==1)/fs],ylim,'r')
    xlim([0 max(data_times)])
    %xlabel('time (s)')
    xticklabels([])
    ylabel('50 Hz relative power')

    nexttile
    plot(steps/fs,ll(:,j))
    hold on
    %plot([steps(ad_on(:,j)==1)/fs,steps(ad_on(:,j)==1)/fs],ylim,'g')
    arrayfun(@(x) xline(x, 'Color', 'g', 'LineWidth', 2), steps(ad_on(:,j)==1)/fs);
    xlim([0 max(data_times)])
    %xlabel('time (s)')
    xticklabels([])
    ylabel('Line length')
   % title(t,sprintf('%s start %1.1f',chs_of_interest{j},start_time))
end

end

%% Start looking for stim and ADs

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

        % 50 Hz filter
        [b, a] = butter(2, [48 52]/(fs/2), 'bandpass');
        filtered_50 = filter(b, a, data);

        % Loop over stim channels
        for j = 1:2
            ch = strcmp(chLabels,stim_chs{j});
            if sum(ch) == 0,error('cannot find stim ch'); end

            data_ch = data(:,ch);
            filtered_50_ch = filtered_50(:,ch);
            rel_50 = sum(filtered_50_ch.^2)/sum(data_ch.^2);

            if rel_50 > stim_start_thresh
                stim_begun(j) = i; % say stim has begun on that channel
            end
        end

    % If stim has begun on both channels and hasn't ended on both channels
    elseif all(stim_begun~=0) && any(stim_ended==0)
        % Look for end of stim
        % 50 Hz filter
        [b, a] = butter(2, [48 52]/(fs/2), 'bandpass');
        filtered_50 = filter(b, a, data);

        % Loop over stim channels
        for j = 1:2
            ch = strcmp(chLabels,stim_chs{j});
            if sum(ch) == 0,error('cannot find stim ch'); end

            data_ch = data(:,ch);
            filtered_50_ch = filtered_50(:,ch);
            rel_50 = sum(filtered_50_ch.^2)/sum(data_ch.^2);

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
                ad_begun = i;
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
            ad_ended = i;

            % Plot stuff
            plot_stims(chLabels,chs_of_interest,data_times,data_stream,...
    stim_begun,stim_ended,ad_begun,ad_ended,ad_ch)
            error('look')

        end

    elseif all(stim_begun~=0) && all(stim_ended~=0) && ad_begun ~=0 && ad_ended ~=0

        % everything done!
    else
        error('what')

    end

end







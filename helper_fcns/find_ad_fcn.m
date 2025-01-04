function T = find_ad_fcn(file_name,start_time,end_time)

%% Parameters

% Display option parameters
do_plots = 1;

% Time chunk parameters
chunkDuration = 0.02; % how long to pull (s)
updateInterval = 0.02; % how much to advance after each pull (s)

% Stim detection parameters
decay = 0.3; % decay rate in adding old signal to new one (looks for repeated signal at 50 Hz)
stimPowerBoost = 1e9; % How big does buffer power need to be
stimOffPower =  1e7; % if drops below this, no longer stim
secs_thresh_stim = 0.3; % How long should it exceed this power
num_thresh_stim = ceil(secs_thresh_stim/chunkDuration);
perc_above_thresh_stim = 0.5; % What percentage should exceed this power?

% Ad detection parameters
n_baseline_all = 100; % How many baselines to take
n_baseline_keep = 50; % How many to keep (only keep first 50 because will assume weirdness right before official stim detection)
stopLookingADSecs = 5; % Stop looking for ADs this long after stim offset 
ad_thresh = 30; % relative power above baseline threshold increased 20->30 1/3
ad_too_high_thresh = 1e4; % if relative power above this, assume artifact
coolDownLastSat = 2; % if ch saturated within this time period, don't look! % reduced 2->1 then back to 2
secs_thresh = 2; % How long does it have the opportunity to get the num above thresh % increased 1->2
num_above_thresh = 10; % How many chunks need to be above threshold to trigger detection 
hfband = [400 500]; % look for high frequency power as an artifact detector 
hfthresh = 1e4; % if hf energy above this, dont count it as being above threshold for AD detection 

% Bad channel parameters
bad_ch_amp = 1e4; % add a bad count if exceeds this outside of stim
n_bad = 10; % if more than this number above bad_ch_amp outside of stim, call it bad
n_bad_reset = 100; % reduce bad ch count by 1 after this many loops
n_reduce_bad = 5;

%% File locs and set path
locations = seizure_termination_paths; % get paths

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% ieeg stuff
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

%% Pull data once for the purpose of getting channel info, etc.
data = download_ieeg_data(file_name,login_name,pwfile,[start_time start_time+1],1);
samplingRate = data.fs;

% channel labels
chLabels = data.chLabels(:,1);
chLabels = decompose_labels(chLabels);
numChannels = length(chLabels);
exclude = find_exclude_chs(chLabels);

% start bad channel counter
bad_ch_counter = zeros(1,numChannels);
n_bad_loop_counter = 0;

aT = data.aT; % this is just for validation (will remove for real time processing)


%% Convert parameters to time samples
chunkSize = round(chunkDuration * samplingRate);
updateSize = round(updateInterval * samplingRate); 
num_thresh = round(secs_thresh/chunkDuration);

%% Initialize variables
% initialize tracker of ad and stim
T = table('Size',[0 4],'VariableTypes',{'cell','cell','double','double'},'VariableNames',{'Type','Channels','OnTime','OffTime'});
curr_row = 0; % tracks AD and detections

% initialize stim detection info
stim_on = 0;
look_for_stim = 1; % are we looking for stim
last_stim_chs = [];
last_stim_off = 1e30;
look_for_sat = 0;

buffer = zeros(chunkSize,numChannels);
last_ones = zeros(num_thresh,numChannels); % counting how many of the last num_thresh werea above threshold
last_ones_stim = zeros(num_thresh_stim,numChannels); % counting how many of the last num_thresh werea above threshold
time_since_last_sat = nan(numChannels,1);
baselines_all = cell(n_baseline_all,1);
for in = 1:n_baseline_all
    baselines_all{in} = zeros(chunkSize,numChannels);
end
baselines = baselines_all(1:n_baseline_keep);
zi = zeros(1,numChannels);

% Look over chunks
%orig_end_time = end_time;
%while orig_end_time > end_time

%% Get ieeg data
data = download_ieeg_data(file_name,login_name,pwfile,[start_time end_time],0); % 1 means get lots of data
values = data.values;

%% High pass filter
% needed for good performance (removes low frequency decay after stim
% artifact). NEED TO IMPLEMENT WITHIN 20 MS LOOP
%[values,zf] = stevefilter(values,zi);
%zi = zf;

% initialize variables for troubleshooting, not needed for final processing
numSamples = size(values,1);
out_samples = 1:updateSize:(numSamples - chunkSize);
all_all_power = nan(length(out_samples),numChannels);
all_all_power_avg = nan(length(out_samples),numChannels);
all_hfp = nan(length(out_samples),numChannels);


count = 0; % troubleshooting thing

% Loop over 20 ms intervals
for startIdx = 1:updateSize:(numSamples - chunkSize)

    % increase the n_bad_loop counteer
    n_bad_loop_counter = n_bad_loop_counter + 1;

    % if you've gone through n_bad_reset loops
    if n_bad_loop_counter == n_bad_reset

        % reduce bad ch counter by one at the reset (min 0)
        bad_ch_counter = max([zeros(1,numChannels);bad_ch_counter-n_reduce_bad*ones(1,numChannels)],[],1);
        n_bad_loop_counter = 0; % and reset the loop counter
    end
    
    
    % Initially, say we're not looking for offset of stim
    look_for_off = 0;
    count = count + 1; % troubleshooting variable


    % define the data chunk
    endIdx = startIdx + chunkSize - 1;
    dataChunk = values(startIdx:endIdx, :);

    % exception handling for ieeg nan periods
    if sum(~isnan(dataChunk),'all') == 0
        dataChunk = zeros(size(dataChunk));
    end

    % replace nans with mean across remaining times
    %
    for ich = 1:numChannels
        dataChunk(isnan(dataChunk(:,ich)),ich) = mean(dataChunk(:,ich),"omitnan");
    end
    %}


    % demean the chunk
    dataChunk = dataChunk - mean(dataChunk,1,"omitnan");
    
    % HPF
    oldDataChunk = dataChunk;
    [dataChunk,zf] = stevefilter(dataChunk,zi);
    zi = zf;

    % remove excluded channels
    dataChunk(exclude) = nan;

    % CAR
    %dataChunk = dataChunk - median(dataChunk,2,"omitnan");

    % Bipolar ref
    %
    [~,bipolarIndices,altBipolarIndices] = find_bipolar_pairs(chLabels,1:length(chLabels));
    newDataChunk = dataChunk;
    newDataChunk(:,~isnan(altBipolarIndices)) = dataChunk(:,~isnan(altBipolarIndices)) - dataChunk(:,altBipolarIndices(~isnan(altBipolarIndices)));
    newDataChunk(:,isnan(altBipolarIndices)) = nan;
    newDataChunk = newDataChunk - mean(newDataChunk,1,"omitnan");
    %{
    newDataChunk = dataChunk;
    for k = 1:length(chLabels)
        r = k == bipolarIndices(:,1);
        if sum(r) == 0, continue; end
        newDataChunk(:,k) = dataChunk(:,bipolarIndices(r,1)) - dataChunk(:,bipolarIndices(r,2));
    end
    %dataChunk = newDataChunk;
    %}

    % If stim is not on, decide whether to look for stim and AD
    if stim_on == 0
        look_for_stim = 1; % look for stim if not current stim
        look_for_off = 0; % don't look for offset if not current stim
        

        
        % Look for AD if within a certain time period of last stim
        if isempty(last_stim_chs)
            look_for_ad = 0;
        %elseif endIdx > last_stim_off + stopLookingAD  
        elseif endIdx/samplingRate + start_time > last_stim_off + stopLookingADSecs
           
            look_for_ad = 0; % don't look for AD if too far after offset of last stim
            look_for_sat = 0;
        else
            look_for_ad = 1; % otherwise, if stim off, look for AD
        end
        %}

    end



    %{
    if ~isempty(last_stim_chs)
    if endIdx/samplingRate + start_time > last_stim_on + coolDownLastOn
        error('what')
    end
    end
    %}

    if stim_on == 1 % if stim is on, look for stim offset and not stim onset
        look_for_off = 1;
        look_for_stim = 0;
        look_for_sat = 1;
    end

    %% Update buffer (used to look for stim)
    buffer(isnan(buffer)) = dataChunk(isnan(buffer)); % exception handling for nans in ieeg
    buffer = buffer*decay+dataChunk; % add dataChunk to last buffer (good for repeated signals)
    buffer_power = sum(buffer.^2,1); % get power
    all_all_power(count,:) = buffer_power;

    %% Look for stim
    if look_for_stim == 1

        % not during stim, look for channels to add to bad ch counter
        bad_ch_counter = bad_ch_counter + sum(newDataChunk > bad_ch_amp,1);
    
        % Get channels meeting amplitude criteria
        chs_above_thresh = buffer_power > stimPowerBoost;

        % Update the tracker of how many channels met this amplitude
        % criteria
        last_ones_stim(1:end-1,:) = last_ones_stim(2:end,:);
        last_ones_stim(end,:) = chs_above_thresh;
    
        % Say stim was detected across a channel if enough in the tracker
        % met the threshold
        detected_stim = sum(last_ones_stim==1,1) > size(last_ones_stim,1)*perc_above_thresh_stim;

        % See how many of these are bipolar channels
        [~,bipolar_indices]= find_bipolar_pairs(chLabels(detected_stim),find(detected_stim));
        
        % if bipolar indices above threshold, then stim is detected
        if ~isempty(bipolar_indices) % stim detected!

            % if more than one bipolar pair above threshold, find the highest
            % amplitude, and say that's the one with stim
            if size(bipolar_indices,1) > 1 
                mean_rel_powers = nan(size(bipolar_indices,1),1);
                for j = 1:size(bipolar_indices,1)
                    mean_rel_powers(j) = mean([buffer_power(bipolar_indices(j,1)),buffer_power(bipolar_indices(j,2))]);
                end
                [~,highest] = max(mean_rel_powers);
                keep_bipolar_indices = bipolar_indices(highest,:);
            else
                keep_bipolar_indices = bipolar_indices;
        
            end
        
            % Update the table tracking stim and AD detections. Each member
            % of the bipolar pair gets a detection
            curr_row = curr_row + 1;
            T.Type{curr_row} = 'stim';
            T.Channels{curr_row} = chLabels{keep_bipolar_indices(1)};
            T.OnTime(curr_row) = endIdx/samplingRate + start_time;
            T.OffTime(curr_row) = nan;
            
            curr_row = curr_row + 1;
            T.Type{curr_row} = 'stim';
            T.Channels{curr_row} = chLabels{keep_bipolar_indices(2)};
            T.OnTime(curr_row) = endIdx/samplingRate + start_time;
            T.OffTime(curr_row) = nan;

            % Reset the tracker of how many met threshold for stim
            % detection
            last_ones_stim(1:end) = 0;

            % Reset the rows of the stim tracker table
            last_stim_rows = [curr_row-1, curr_row];
            stim_on = 1; % say stim is on
            last_stim_on = T.OnTime(curr_row);
            look_for_ad = 0; % say we should not look for ADs
            ad_chs_this_stim = []; % reset what channels had ADs for this stim period
            last_stim_chs = [keep_bipolar_indices(1) keep_bipolar_indices(2)]; % set these channels to be the last stim channels


            % Get multiple baselines (for comparison for AD detection)
            %{
            baselines = cell(bl_chunks(2),1);
            for j = 1:bl_chunks(2)
                curr_bl = values(startIdx - chunkSize*bl_chunks(1) + (j-1)*chunkSize + 1: ...
                    startIdx - chunkSize*bl_chunks(1) + (j)*chunkSize,:);
                curr_bl = curr_bl - nanmean(curr_bl,1);
                baselines{j} = curr_bl;
            end
            %}

            % Set baseline to only be the first n_baseline_keep (given that
            % more recent ones will be close to threshold)
            baselines = baselines_all(1:n_baseline_keep);


        end

        baselines_all(1:end-1) = baselines_all(2:end);
        %baselines_all(end) = {dataChunk};
        baselines_all(end) = {newDataChunk};
    
        
    end

    % if looking for stim offset
    if look_for_off
        stim_chs = [keep_bipolar_indices(1) keep_bipolar_indices(2)];
        buffer_power_stim = buffer_power(stim_chs);
        if mean(buffer_power_stim) < stimOffPower
            stim_on = 0; % stim is off
            T.OffTime(last_stim_rows(1)) = endIdx/samplingRate + start_time;
            T.OffTime(last_stim_rows(2)) = endIdx/samplingRate + start_time;
            %last_stim_off = endIdx;
            last_stim_off = T.OffTime(last_stim_rows(1));
        end
    end

    %% Look for ad
    
        %error('what')
    if look_for_sat || look_for_ad
        % calc power
        %power = measure_power(dataChunk,freq_range,samplingRate,chLabels);
        power = measure_power(newDataChunk);
        power_avg = power;

        

        % Get the power for the baselines using the same approach
        baseline_power = cellfun(@(x) measure_power(x), baselines, ...
            'UniformOutput',false);
        baseline_power = cell2mat(baseline_power);
        baseline_power = mean(baseline_power,1); % Take the mean across the baseline
        baseline_power_avg = baseline_power;


        % get the relative power compared to mean of baseline
        rel_power_avg = power_avg./baseline_power_avg;

        % have a special check for likely artifact        
        above_thresh = rel_power_avg > ad_thresh & rel_power_avg < ad_too_high_thresh;

        % Look for high frequency power - make zero if too high
        high_freq_power = measure_power(newDataChunk,hfband,samplingRate);
        all_hfp(count,:) = high_freq_power;
        too_much_high_freq = high_freq_power > hfthresh;
        above_thresh(too_much_high_freq) = 0;

        % saturated?
        saturated = rel_power_avg > ad_too_high_thresh;
        time_since_last_sat(saturated) = endIdx/samplingRate + start_time; % say it saturated now
        
        % reset above thresh designation for those within the saturated
        % cool down
        above_thresh(endIdx/samplingRate + start_time-time_since_last_sat<coolDownLastSat) = 0;

        % For troubleshooting
        all_all_power(count,:) = power;
        all_all_power_avg(count,:) = rel_power_avg;
    end

    if look_for_ad 
        % Reset tracker for how many were above thresh
        last_ones(1:end-1,:) = last_ones(2:end,:);
        last_ones(end,:) = above_thresh;


        % Decide if enough above thresh
        %detected_ad = sum(last_ones==1,1) > size(last_ones,1)*perc_above_thresh;
        detected_ad = sum(last_ones==1,1) > num_above_thresh;

        

        % make it zero if it's not an ad look channel (anything not on the stim electrode, excluding stim contacts!)
        chs_to_look = ad_look_chs(chLabels,altBipolarIndices,last_stim_chs,file_name);
        detected_ad(chs_to_look == 0) = 0;

        % make zero if bad ch counter above limit
        detected_ad(bad_ch_counter > n_bad) = 0;

        % if detection, add it
        if any(detected_ad,'all')      
            %error('what')

            % find the detected ones
            detected_ones = find(detected_ad);
            for k = 1:length(detected_ones)
                % skip if already an AD for this channel and this stim
                if ismember(detected_ones(k),ad_chs_this_stim)
                    continue
                end

                % add it to the ad table
                curr_row = curr_row + 1;
                T.Type{curr_row} = 'AD';
                T.OnTime(curr_row) = endIdx/samplingRate + start_time;
                T.Channels{curr_row} = chLabels{detected_ones(k)};
                T.OffTime(curr_row) = nan;
                last_ad_row = curr_row;
                ad_chs_this_stim = [ad_chs_this_stim;detected_ones(k)];
            end

            
            last_ones(:,ad_chs_this_stim) = 0; % reset all to zero

        end


    else
        % reset AD thresh counter if not looking for ADs
        last_ones = zeros(num_thresh,numChannels);

    end


         
end

% Plot the EEG data with markers for stim and afterdischarge detections (as before)

if do_plots
    fig_rows = (T.OnTime > start_time & (isnan(T.OffTime) | T.OffTime < end_time));
if sum(fig_rows) == 0
    fprintf('\nNo stim detections or ad detections for %1.1f\n',start_time);
else

    % get bipolar
    [~,bipolarIndices,altBipolarIndices] = find_bipolar_pairs(chLabels,1:length(chLabels));
    biValues = values;
    biValues(:,~isnan(altBipolarIndices)) = values(:,~isnan(altBipolarIndices)) - values(:,altBipolarIndices(~isnan(altBipolarIndices)));
    biValues(:,isnan(altBipolarIndices)) = nan;
    biValues = biValues - mean(biValues,1,"omitnan");

    N = sum(ismember(chLabels,T.Channels(fig_rows)));
    rows = ceil(sqrt(N));
    cols = ceil(N / rows);
    
    
    figure;
    t = tiledlayout(rows, cols);
    for ch = 1:numChannels
    
        if ismember(chLabels(ch),T.Channels(fig_rows))
            nexttile
            % Plot the raw EEG data for the current channel
            plot(linspace(start_time,end_time,size(values,1)),biValues(:, ch), 'k'); % Plot in black for baseline EEG
            hold on

            % get rows of stim detections for that channel
            stim_rows = find(strcmp(chLabels(ch),T.Channels) & strcmp(T.Type,'stim') & fig_rows);
            stim_on_times = T.OnTime(stim_rows);
            stim_off_times = T.OffTime(stim_rows);

            if ~isempty(stim_rows)
                stim_des = 'S';
            else
                stim_des = '';
            end
            
            % Overlay stimulation detections in red
            for j = 1:length(stim_on_times)
                plot([stim_on_times(j) stim_on_times(j)],ylim,'b')
                plot([stim_off_times(j) stim_off_times(j)],ylim,'r')
            end

            % (Troubleshooting) plot the baseline times
            yl = ylim;
            for j = 1:length(stim_on_times)
                plot([stim_on_times(j) - chunkDuration*n_baseline_all stim_on_times(j) - chunkDuration*n_baseline_all+chunkDuration*n_baseline_keep],...
                    [yl(1)+0.9*(yl(2)-yl(1)) yl(1)+0.9*(yl(2)-yl(1))],'r')
            end
            
            % Overlay the closed relay annotations
            rows_in_range = find(aT.Start > start_time & aT.Start < end_time & contains(aT.Type,'Closed relay'));
            for i = 1:length(rows_in_range)
                r = rows_in_range(i);
                new_time = aT.Start(r);
                plot([new_time new_time],ylim,'m')
                
                matches = regexp(aT.Type{r}, 'Closed relay to ([A-Z]+\d+) and ([A-Z]+\d+)', 'tokens');
                if ~isempty(matches)
                    firstString = matches{1}{1};  % First captured string
                    secondString = matches{1}{2}; % Second captured string
                    yl = get(gca,'ylim');
                    text(new_time,yl(2),sprintf('%s-%s',firstString,secondString),'Color','m')
                end
                
    
    
            end
    
            % Overlay ad detections
            % get rows of ad detections for that channel
            ad_rows = find(strcmp(chLabels(ch),T.Channels) & strcmp(T.Type,'AD')& fig_rows);
            if ~isempty(ad_rows)
                ad_times = T.OnTime(ad_rows);
                for j = 1:length(ad_times)
                    plot([ad_times(j) ad_times(j)],ylim,'g')
                end

                
                ad_des = 'AD';
                
            else
                ad_des = '';
            end
    
    
            % Add titles and labels for clarity
            title(sprintf('%s %s %s',chLabels{ch},stim_des,ad_des));
            xlabel('Time (s)');
            ax = gca; ax.XRuler.Exponent = 0;
            ylabel('Amplitude');
            
          
        
        end
    end
    title(t,sprintf('%1.1f s',start_time))
end
%pause
end

if 0
curr_labs = {'LE6'};
for i = 1:length(curr_labs)
    
    ch = strcmp(chLabels,curr_labs{i}); %RI3
    if sum(~isnan(values(:,ch))) == 0, continue; end
    figure
    nexttile
    plot(linspace(start_time,end_time,size(values,1)),values(:,ch)-median(values(:,~exclude),2,"omitnan"))
    nexttile
    alpha = 0.99;
    [xf,zf] = filter(1-alpha,[1, -alpha],values,zeros(1,numChannels));
    xout = values - xf;
    newValues = xout;
    newValues(:,~isnan(altBipolarIndices)) = newValues(:,~isnan(altBipolarIndices)) - newValues(:,altBipolarIndices(~isnan(altBipolarIndices)));
    newValues = newValues - mean(newValues,1,"omitnan");
    plot(linspace(start_time,end_time,size(newValues,1)),newValues(:,ch)-median(newValues(:,~exclude),2,"omitnan"))
    nexttile
    plot(linspace(start_time,end_time,size(all_all_power_avg,1)),all_all_power_avg(:,ch))
    nexttile
    spectrogram(newValues(:,ch),100,80,100,samplingRate,'yaxis')
    nexttile
    plot(linspace(start_time,end_time,size(all_hfp,1)),all_hfp(:,ch))
end
end





%T
%return
%fprintf('\nMove on?\n');
%pause
%fprintf('\nMoving on\n');

%start_time = start_time + duration;
%end_time = end_time + duration;
%end

% clean it up
%ad_empty = cellfun(@(x) isempty(x),T.AD);
%T.AD_display = T.AD;
%T.AD_display(ad_empty) = {{''}};
%T.AD_display = cellfun(@(x) strjoin(x,','),T.AD_display,'UniformOutput',false);

end

%figure
%imagesc(all_all_rel_power')

function chs_to_look = ad_look_chs(chLabels,altBipolarIndices,stimChs,file_name)

%% Grab all channels
all_chs = (1:length(chLabels))';
chs_to_look = logical(zeros(length(chLabels),1));

%% Find channels on same electrode as stim channel
is_stim_ch = (ismember(all_chs,stimChs) | ismember(altBipolarIndices,stimChs));
% get the stim channel labels
stimChLabels = chLabels(stimChs);

% Take the first one
stimChLabel = stimChLabels{1};

% Get the letter part
match = regexp(stimChLabel, '([A-Z]+)(\d+)', 'tokens');
letterPart = match{1}{1};

% Find all channels that start with the same letters
checkLabel = @(x) ...
    (length(x) >= length(letterPart)) && strcmp(letterPart, x(1:length(letterPart)));

same_elec = cellfun(checkLabel, chLabels);

%% Find channels on contact <11 (to avoid out of brain)
% Apply regexp to each element of the cell array to find numbers at the end
% Use regexp to find a trailing sequence of digits (if any) for each string
numericParts = cellfun(@(s) regexp(s, '\d+$', 'match'), chLabels, 'UniformOutput', false);

% Initialize a numeric array filled with NaN
numericArray = nan(size(chLabels));

% Find which elements are not empty
notEmptyIdx = ~cellfun(@isempty, numericParts);

% For those with matches, convert the matched string to a double
numericArray(notEmptyIdx) = cellfun(@(x) str2double(x{1}), numericParts(notEmptyIdx));


less_than_eleven = numericArray < 11; % higher numbered outside brain and susceptible to noise;

%% Define what to look for
% Look if it's on same elec as stim channels but is not itself a stim ch
if strcmp(file_name,'HUP260_phaseII')
    chs_to_look(same_elec & ~is_stim_ch) = 1;
else
    chs_to_look(same_elec & ~is_stim_ch & less_than_eleven) = 1;
end
%chs_to_look(same_elec) = 1;
%chs_to_look(~is_stim_ch) = 1;

% exclude excludable channels
exc = find_exclude_chs(chLabels);
chs_to_look(exc) = 0;

end


function exclude = find_exclude_chs(chLabels)

excluded = {'C3','C4','CZ','EKG1','EKG2','FZ','LOC','ROC'};
exclude = ismember(chLabels,excluded);

end


function [bipolarPairs,bipolarIndices,altBipolarIndices] = find_bipolar_pairs(contacts,contactIndices)

% Initialize cell array to store bipolar pairs
bipolarPairs = {};
bipolarIndices = [];
altBipolarIndices = nan(length(contacts),1);

% Loop through each contact to extract letter and number parts
for i = 1:numel(contacts)
    % Parse the current contact
    match = regexp(contacts{i}, '([A-Z]+)(\d+)', 'tokens');
    if ~isempty(match)
        letterPart = match{1}{1};               % Extract letter part
        numberPart = str2double(match{1}{2});    % Extract number part as double

        % Check for the next contact with the same letter part and incremented number
        nextContact = sprintf('%s%d', letterPart, numberPart + 1);
        if any(strcmp(contacts, nextContact))
            % If next contact exists, add the pair to bipolarPairs
            bipolarPairs = [bipolarPairs; {contacts{i}, nextContact}];
            bipolarIndices = [bipolarIndices; contactIndices(i) contactIndices(find(strcmp(contacts,nextContact)))];
            altBipolarIndices(i) = find(strcmp(contacts,nextContact));
        end
    end
end


end



function p = measure_power(data,freq_band,fs)

meanValues = mean(data,1,'omitnan');
nanIndices = isnan(data);
for col = 1:size(data, 2)
    data(nanIndices(:, col), col) = meanValues(col); % this takes a LOONG time. why am I doing this????
end

if exist('freq_band','var') == 0

    % just do abs power
    p = sum(data.^2,1);
else
    data(isnan(data)) = 0;
    p = bandpower(data,fs,freq_band);
end

end

function [xout,zf] = stevefilter(xin,zi)
    alpha = 0.99; 
    [xf,zf] = filter(1-alpha,[1, -alpha],xin,zi);
    xout = xin - xf;
end


%{
to do
- for some reason if I look at a long period of time (100 s) starting at
22368 it falsely detects ADs, but not a shorter time. Must have sometihng
to do with resetting with stim

Good times to check:
- HUP218_HFS: 20793.72 look

- really need to find any stim and start looking as soon as stim is done

%}

clear

% Parameters
do_plots = 1;
chunkDuration = 0.02; 
updateInterval = 0.02; 
coolDownPeriodSecsAD = 5;
stopLookingADSecs = 10;
decay = 0.3;
stimPowerBoost = 1e9;
freq_range = [5 40];
ad_thresh = 50;
secs_thresh = 1;
secs_thresh_stim = 0.5;
num_thresh = round(secs_thresh/chunkDuration);
num_thresh_stim = ceil(secs_thresh_stim/chunkDuration);
bl_chunks = [20 10];
perc_above_thresh = 0.2;
perc_above_thresh_stim = 0.5;

%% Main
file_name = 'HUP225_HFS';%'HUP218_HFS';
start_time = 13697.85;%14057.47;%20718.06;
duration = 300;
end_time = start_time + duration;


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

% initialize tracker of ad and stim
T = table('Size',[0 4],'VariableTypes',{'cell','cell','double','double'},'VariableNames',{'Type','Channels','OnTime','OffTime'});
curr_row = 0;

stim_on = 0;
look_for_stim = 1;
last_stim_chs = [];

while 1
%% Get ieeg data
data = download_ieeg_data(file_name,login_name,pwfile,[start_time end_time],1); % 1 means get lots of data
values = data.values;
times = linspace(0,duration,size(values,1));
chLabels = data.chLabels(:,1);
chLabels = decompose_labels(chLabels);
samplingRate = data.fs;
[numSamples, numChannels] = size(values);
aT = data.aT;

exclude = find_exclude_chs(chLabels);


% Convert parameters to time samples
chunkSize = round(chunkDuration * samplingRate);
updateSize = round(updateInterval * samplingRate); 
cooldownPeriodAD = round(coolDownPeriodSecsAD * samplingRate); 
stopLookingAD = round(stopLookingADSecs*samplingRate);

out_samples = 1:updateSize:(numSamples - chunkSize);





buffer = zeros(chunkSize,numChannels);
all_all_power = nan(length(out_samples),numChannels);
all_all_power_avg = nan(length(out_samples),numChannels);
last_ones = zeros(num_thresh,numChannels); % counting how many of the last num_thresh werea above threshold
last_ones_stim = zeros(num_thresh_stim,numChannels); % counting how many of the last num_thresh werea above threshold


count = 0; % troubleshooting thing

for startIdx = 1:updateSize:(numSamples - chunkSize)
    
    look_for_off = 0;
    

    count = count + 1;

    endIdx = startIdx + chunkSize - 1;
    dataChunk = values(startIdx:endIdx, :);

    % demean 
    old_dataChunk = dataChunk;
    dataChunk = dataChunk - mean(dataChunk,1);

    % remove excluded channels
    dataChunk(exclude) = nan;



    if stim_on == 0
        look_for_stim = 1;
        look_for_off = 0;
        % Look for AD if within a certain time period of last stim
        if isempty(last_stim_chs)
            look_for_ad = 0;
        elseif endIdx > last_stim_off + stopLookingAD  
            look_for_ad = 0;
         
        else
            look_for_ad = 1;
        end

    end

    if stim_on == 1
        look_for_off = 1;
        look_for_stim = 0;
    end

    % Update buffer
    buffer = buffer*decay + dataChunk;
    buffer_power = sum(buffer.^2,1);
    all_all_power(count,:) = buffer_power;

    %% Look for stim
    if look_for_stim == 1
    
        % Get channels meeting amplitude criteria
        chs_above_thresh = buffer_power > stimPowerBoost;
        last_ones_stim(1:end-1,:) = last_ones_stim(2:end,:);
        last_ones_stim(end,:) = chs_above_thresh;
    
        detected_stim = sum(last_ones_stim==1,1) > size(last_ones_stim,1)*perc_above_thresh_stim;

        % See how many of these are bipolar channels
        [~,bipolar_indices]= find_bipolar_pairs(chLabels(detected_stim),find(detected_stim));
        %[~,bipolar_indices]= find_bipolar_pairs(chLabels(chs_above_thresh),find(chs_above_thresh));
        
        % skip if no bipolar indices above threshold
        if ~isempty(bipolar_indices)

            % if more than one bipolar pair above threshold, find the highest
            % amplitude
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
        
            curr_row = curr_row + 1;
            T.Type{curr_row} = 'stim';
            T.Channels{curr_row} = chLabels{keep_bipolar_indices(1)};
            %T.Channels{curr_row} = sprintf('%s-%s',chLabels{keep_bipolar_indices(1)},chLabels{keep_bipolar_indices(2)});
            T.OnTime(curr_row) = endIdx/samplingRate + start_time;
            T.OffTime(curr_row) = nan;
            
            curr_row = curr_row + 1;
            T.Type{curr_row} = 'stim';
            T.Channels{curr_row} = chLabels{keep_bipolar_indices(2)};
            %T.Channels{curr_row} = sprintf('%s-%s',chLabels{keep_bipolar_indices(1)},chLabels{keep_bipolar_indices(2)});
            T.OnTime(curr_row) = endIdx/samplingRate + start_time;
            T.OffTime(curr_row) = nan;

            last_ones_stim(1:end) = 0;

            last_stim_rows = [curr_row-1, curr_row];
            stim_on = 1;
            look_for_ad = 0;
            ad_chs_this_stim = [];
            last_stim_chs = [keep_bipolar_indices(1) keep_bipolar_indices(2)];


            % Get multiple baseline
            baselines = cell(bl_chunks(2),1);
            for j = 1:bl_chunks(2)
                curr_bl = values(startIdx - chunkSize*bl_chunks(1) + (j-1)*chunkSize + 1: ...
                    startIdx - chunkSize*bl_chunks(1) + (j)*chunkSize,:);
                curr_bl = curr_bl - nanmean(curr_bl,1);
                baselines{j} = curr_bl;
            end

        end
    
        
    end

    if look_for_off
        stim_chs = [keep_bipolar_indices(1) keep_bipolar_indices(2)];
        buffer_power_stim = buffer_power(stim_chs);
        if mean(buffer_power_stim) < 1e7
            stim_on = 0;
            T.OffTime(last_stim_rows(1)) = endIdx/samplingRate + start_time;
            T.OffTime(last_stim_rows(2)) = endIdx/samplingRate + start_time;
            last_stim_off = endIdx;
        end
    end

    %% Look for ad
    if look_for_ad % issue here
        
        % calc power
        power = measure_power(dataChunk,freq_range,samplingRate);
        power_avg = power;

        baseline_power = cellfun(@(x) measure_power(x,freq_range,samplingRate), baselines, ...
            'UniformOutput',false);
        baseline_power = cell2mat(baseline_power);
        baseline_power = mean(baseline_power,1); % Take the mean across the baseline
        baseline_power_avg = baseline_power;



        rel_power_avg = power_avg./baseline_power_avg;

        % have a special check for likely artifact
        
        above_thresh = rel_power_avg > ad_thresh & rel_power_avg < 1e4;
        last_ones(1:end-1,:) = last_ones(2:end,:);
        last_ones(end,:) = above_thresh;
        %{
        if endIdx/samplingRate > 10.5
            endIdx/samplingRate
            rel_power_avg(177)
            last_ones(:,177)
            pause
        end
        %}

        all_all_power(count,:) = power;
        all_all_power_avg(count,:) = rel_power_avg;

        % Decide if enough above thresh
        %detected_ad = num_above_thresh > num_thresh;
        % majority above thresh
        detected_ad = sum(last_ones==1,1) > size(last_ones,1)*perc_above_thresh;

        % make it zero if it's not an ad look channel
        chs_to_look = ad_look_chs(chLabels,last_stim_chs);
        detected_ad(chs_to_look == 0) = 0;

        % if detection, add it
        if any(detected_ad,'all')      
            %error('what')
            last_ones(:) = 0;

            detected_ones = find(detected_ad);
            for k = 1:length(detected_ones)
                % ADD SOMETHING TO SKIP IF ALREADY AN AD FOR THIS STIM
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

        end


    else
        % reset AD thresh counter
        last_ones = zeros(num_thresh,numChannels);

    end


         
end

% Plot the EEG data with markers for stim and afterdischarge detections (as before)

if do_plots
    fig_rows = (T.OnTime > start_time & (isnan(T.OffTime) | T.OffTime < end_time));
if sum(fig_rows) == 0
    fprintf('\nNo stim detections or ad detections for %1.1f\n',start_time);
else

    N = sum(ismember(chLabels,T.Channels(fig_rows)));
    rows = ceil(sqrt(N));
    cols = ceil(N / rows);
    
    
    figure;
    t = tiledlayout(rows, cols);
    for ch = 1:numChannels
    
        if ismember(chLabels(ch),T.Channels(fig_rows))
            nexttile
            % Plot the raw EEG data for the current channel
            plot(linspace(start_time,end_time,size(values,1)),values(:, ch), 'k'); % Plot in black for baseline EEG
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
curr_labs = {'RR4'};
for i = 1:length(curr_labs)
    ch = strcmp(chLabels,curr_labs{i}); %RI3
    figure
    nexttile
    plot(linspace(0,duration,size(values,1)),values(:,ch))
    nexttile
    plot(linspace(0,duration,size(all_all_power_avg,1)),all_all_power_avg(:,ch))
    nexttile
    spectrogram(values(:,ch),100,80,100,samplingRate,'yaxis')
end
end





T
%return
fprintf('\nMove on?\n');
pause
fprintf('\nMoving on\n');

start_time = start_time + duration;
end_time = end_time + duration;
end

% clean it up
ad_empty = cellfun(@(x) isempty(x),T.AD);
T.AD_display = T.AD;
T.AD_display(ad_empty) = {{''}};
T.AD_display = cellfun(@(x) strjoin(x,','),T.AD_display,'UniformOutput',false);


%figure
%imagesc(all_all_rel_power')

function chs_to_look = ad_look_chs(chLabels,stimChs)

all_chs = 1:length(chLabels);
is_stim_ch = ismember(all_chs,stimChs)';
chs_to_look = logical(zeros(length(chLabels),1));

% get the stim channel labels
stimChLabels = chLabels(stimChs);

% Take the first one
stimChLabel = stimChLabels{1};

% Get the letter part
match = regexp(stimChLabel, '([A-Z]+)(\d+)', 'tokens');
letterPart = match{1}{1};

% Find all channels that start with the same letters
same_elec = cellfun(@(x) strcmp(letterPart,x(1:length(letterPart))),chLabels);

% Look if it's on same elec as stim channels but is not itself a stim ch
chs_to_look(same_elec & ~is_stim_ch) = 1;
%chs_to_look(same_elec) = 1;

% exclude excludable channels
exc = find_exclude_chs(chLabels);
chs_to_look(exc) = 0;

end


function exclude = find_exclude_chs(chLabels)

excluded = {'C3','C4','CZ','EKG1','EKG2','FZ','LOC','ROC'};
exclude = ismember(chLabels,excluded);

end


function [bipolarPairs,bipolarIndices] = find_bipolar_pairs(contacts,contactIndices)

% Initialize cell array to store bipolar pairs
bipolarPairs = {};
bipolarIndices = [];

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
        end
    end
end


end



function ll = measure_ll(data)

ll = sum(abs(diff(data,1,1)));

end

function p = measure_power(data,freq_range,fs)

meanValues = mean(data,1,'omitnan');
nanIndices = isnan(data);
for col = 1:size(data, 2)
    data(nanIndices(:, col), col) = meanValues(col);
end
p = bandpower(data,fs,freq_range);

% rel power
%p = p./sum(data.^2,1);

end
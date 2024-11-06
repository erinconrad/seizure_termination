%{
Pseudocode:
- Start during a resting period of EEG
- Continuously update:
    - Pull in a time chunk
    - See if stimulation has started, and where (and keep note of this)
    - See if ADs happening, and where (and make a note of this)
    

%}

clear

% Parameters
chunkDuration = 0.02; 
updateInterval = 0.02; 
coolDownPeriodSecs = 6; % 5s
decay = 0.7;
powerBoost = 5e9;%1e3;


%% Main
file_name = 'HUP218_HFS';
start_time = 20793.72; 
duration = 150;
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
cooldownPeriod = round(coolDownPeriodSecs * samplingRate); 
last_detection = -5*samplingRate;

out_samples = 1:updateSize:(numSamples - chunkSize);

detections = [];
buffer = zeros(chunkSize,numChannels);
all_all_power = nan(length(out_samples),numChannels);

count = 0; % troubleshooting thing

for startIdx = 1:updateSize:(numSamples - chunkSize)
    count = count + 1;

    endIdx = startIdx + chunkSize - 1;
    dataChunk = values(startIdx:endIdx, :);

    % demean 
    dataChunk = dataChunk - mean(dataChunk,1);

    % remove excluded channels
    dataChunk(exclude) = nan;

    % Skip if it's in the cool down period
    if endIdx < last_detection + cooldownPeriod
        continue;
    end

    % Update buffer
    buffer = buffer*decay + dataChunk;
    buffer_power = sum(buffer.^2,1);

    all_all_power(count,:) = buffer_power;

    % Get channels meeting amplitude criteria
    chs_above_thresh = buffer_power > powerBoost;

    % See how many of these are bipolar channels
    [bipolar_above_thresh,bipolar_indices]= find_bipolar_pairs(chLabels(chs_above_thresh),find(chs_above_thresh));
    
    % skip if no bipolar indices above threshold
    if isempty(bipolar_indices)
        continue; 
    end

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

    detections = [detections; endIdx keep_bipolar_indices(1); endIdx keep_bipolar_indices(2)];
    last_detection = endIdx;

     
    
end

% Plot the EEG data with markers for stim and afterdischarge detections (as before)
figure;
for ch = 1:numChannels

    if ismember(ch,detections(:,2))
        nexttile
        % Plot the raw EEG data for the current channel
        plot(linspace(0,duration,size(values,1)),values(:, ch), 'k'); % Plot in black for baseline EEG
        hold on
        
        % Overlay stimulation detections in red
        d = ismember(detections(:,2),ch);
        stimIdx = detections(d,1);
        %plot(times(stimIdx), values(stimIdx, ch), 'r.', 'MarkerSize', 10);
        for j = 1:length(stimIdx)
            plot([times(stimIdx(j)) times(stimIdx(j))],ylim,'r')
        end
        
        % Overlay the closed relay annotations
        rows_in_range = find(aT.Start > start_time & aT.Start < end_time & contains(aT.Type,'Closed relay'));
        for i = 1:length(rows_in_range)
            r = rows_in_range(i);
            new_time = aT.Start(r)-start_time;
            plot([new_time new_time],ylim,'b')
            
            matches = regexp(aT.Type{r}, 'Closed relay to ([A-Z]+\d+) and ([A-Z]+\d+)', 'tokens');
            if ~isempty(matches)
                firstString = matches{1}{1};  % First captured string
                secondString = matches{1}{2}; % Second captured string
                yl = get(gca,'ylim');
                text(new_time,yl(2),sprintf('%s-%s',firstString,secondString),'Color','b')
            end
            


        end

        % Add titles and labels for clarity
        title(['Channel ' chLabels(ch)]);
        xlabel('Samples');
        ylabel('Amplitude');
        
      
    
    end
end

pause
start_time = start_time + duration;
end_time = end_time + duration;
end

%figure
%imagesc(all_all_rel_power')



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
function firstBelow = get_precise_offset(file_name,ad_time,ad_ch,stim_times)

%{
The purpose of this function is to take the ieeg filename and time of an
afterdischarge detection and to get the precise offset time
%}

do_plot= 0;

%% Parameters
pre_time = 20;
post_time = 30;
baseline_rel_stim_on = [-2 -1]; % take the one second before stim onset to get baseline
ad_thresh_lower = 5;
chunkDuration = 0.02;
secsLower = 0.5;
percLower = 0.8;

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

%% Pull data
start_time = ad_time-pre_time;
end_time = ad_time + post_time;
data = download_ieeg_data(file_name,login_name,pwfile,[start_time end_time],1);
fs = data.fs;
chunkDurIdx = round(chunkDuration*fs);

chLabels = data.chLabels(:,1);
chLabels = decompose_labels(chLabels);
nchs = length(chLabels);
values = data.values;

%% high pass filter
zi = zeros(1,nchs);
values = stevefilter(values,zi);

%% Get the specific channel
% get bipolar pairs
[~,bipolarIndices,altBipolarIndices] = find_bipolar_pairs(chLabels,1:length(chLabels));
newValues(:,~isnan(altBipolarIndices)) = values(:,~isnan(altBipolarIndices)) - values(:,altBipolarIndices(~isnan(altBipolarIndices)));
newValues(:,isnan(altBipolarIndices)) = nan;
values = newValues;

% restrict to ch
ch_idx = strcmp(chLabels,ad_ch);
values = values(:,ch_idx);
values = values - nanmedian(values);

%% Establish the baseline
baseline_times = [stim_times(1) + baseline_rel_stim_on(1), stim_times(1) + baseline_rel_stim_on(2)];
baseline_times_rel_start = [baseline_times(1) - start_time, baseline_times(2) - start_time];
baseline_idx = round(baseline_times_rel_start(1)*fs):round(baseline_times_rel_start(2)*fs);
values_bl = values(baseline_idx);
n_baseline_chunks = floor(length(values_bl)/chunkDurIdx);
baseline = nan(n_baseline_chunks,1);
for i = 1:n_baseline_chunks
    idx = chunkDurIdx*(i-1)+1:chunkDurIdx*i;
    baseline(i) = sum(values_bl(idx).^2);
end

% get avg power across baseline
avg_bl_power = median(baseline);

%% Get power in 20 ms chunks starting at ad detection time
test_times = [ad_time,end_time];
test_times_rel_start = test_times - start_time;
test_idx = round(test_times_rel_start(1)*fs):round(test_times_rel_start(2)*fs);
values_test = values(test_idx);
n_chunks = floor(length(values_test)/chunkDurIdx);
chunk_times = nan(n_chunks,1);
chunks = nan(n_chunks,1);
for i = 1:n_chunks
    idx = chunkDurIdx*(i-1)+1:chunkDurIdx*i;
    chunks(i) = sum(values_test(idx).^2)/avg_bl_power;
    chunk_times(i) = ad_time + (i-1)*chunkDuration;
end
belowThresh = chunks < ad_thresh_lower; % binary for each chunk if below thresh or not

%% Find point at which drop below thresh for enough percentage of enough time
numBelow = round(secsLower/chunkDuration);
enoughBelow = zeros(n_chunks,1);
for i = numBelow+1:n_chunks
    enoughBelow(i) = sum(belowThresh(i-numBelow:i))/length(belowThresh(i-numBelow:i)) > percLower;
end

if sum(enoughBelow) == 0
    firstBelow = nan;
else
    firstBelowIdx = find(enoughBelow==1);
    firstBelowIdx = firstBelowIdx(1);
    firstBelow = chunk_times(firstBelowIdx);

    % realign to start
    firstBelow = firstBelow - secsLower;
end

if do_plot

    figure
    set(gcf,'Position',[1 1 1400 600])
    tiledlayout(2,1)
    nexttile
    plot(linspace(start_time,end_time,length(values)),values,'k')
    hold on
    plot([ad_time ad_time],get(gca,'ylim'),'g--')
    plot([firstBelow firstBelow],get(gca,'ylim'),'r--')
    plot([baseline_times(1) baseline_times(1)],get(gca,'ylim'),'b--')
    plot([baseline_times(2) baseline_times(2)],get(gca,'ylim'),'b--')
    title(sprintf('%s %s ad times %1.1f %1.1f',file_name,ad_ch,ad_time,firstBelow))

    nexttile
    plot(chunk_times,chunks)
    hold on
    plot([ad_time ad_time],get(gca,'ylim'),'g--')
    plot([firstBelow firstBelow],get(gca,'ylim'),'r--')
    pause

end

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

function [xout,zf] = stevefilter(xin,zi)
    alpha = 0.90;
    [xf,zf] = filter(1-alpha,[1, -alpha],xin,zi);
    xout = xin - xf;
end
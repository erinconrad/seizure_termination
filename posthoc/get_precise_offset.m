function [offsetTime,onsetTime] = get_precise_offset(file_name,ad_time,ad_ch,stim_times)

%{
The purpose of this function is to take the ieeg filename and time of an
afterdischarge detection and to get the precise offset time
%}

do_plot= 0;

%% Parameters
pre_time = 20;
post_time = 30;
baseline_rel_stim_on = [-2 -1]; % take the one second before stim onset to get baseline
chunkDuration = 0.02;
ad_thresh_lower = 5;
secsLower = 0.5;
percLower = 0.8;

ad_thresh_higher = 10;
secsHigher = 0.1;
percHigher = 0.8;



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
    offsetTime = firstBelow;
else
    firstBelowIdx = find(enoughBelow==1);
    firstBelowIdx = firstBelowIdx(1);
    firstBelow = chunk_times(firstBelowIdx);

    % realign to start
    offsetTime = firstBelow - secsLower;
end

%% Also try to get precise onset

% Power in 20 ms chunks starting at stim offset
on_times = [stim_times(2),end_time];
on_times_rel_start = on_times - start_time;
on_idx = round(on_times_rel_start(1)*fs):round(on_times_rel_start(2)*fs);
values_on = values(on_idx);
n_on_chunks = floor(length(values_on)/chunkDurIdx);
on_chunk_times = nan(n_on_chunks,1);
on_chunks = nan(n_on_chunks,1);
for i = 1:n_on_chunks
    idx = chunkDurIdx*(i-1)+1:chunkDurIdx*i;
    on_chunks(i) = sum(values_on(idx).^2)/avg_bl_power;
    on_chunk_times(i) = stim_times(2) + (i-1)*chunkDuration;
end

aboveThresh = on_chunks > ad_thresh_higher;

% Find point at which go above thresh for enough percentage of enough time
numAbove = round(secsHigher/chunkDuration);
enoughAbove = zeros(n_on_chunks,1);
for i = numAbove+1:n_on_chunks
    enoughAbove(i) = sum(aboveThresh(i-numAbove:i))/length(aboveThresh(i-numAbove:i)) > percHigher;
end

if sum(enoughAbove) == 0
    firstAbove = nan;
    onsetTime = firstAbove;
else
    firstAboveIdx = find(enoughAbove==1);
    firstAboveIdx = firstAboveIdx(1);
    firstAbove = on_chunk_times(firstAboveIdx);

    % realign to start
    onsetTime = firstAbove - secsHigher;
end


if do_plot

    figure
    set(gcf,'Position',[1 1 1400 1000])
    tiledlayout(3,1)
    nexttile
    plot(linspace(start_time,end_time,length(values)),values,'k')
    hold on
    plot([ad_time ad_time],get(gca,'ylim'),'g--')
    plot([onsetTime onsetTime],get(gca,'ylim'),'c--')
    plot([offsetTime offsetTime],get(gca,'ylim'),'r--')
    plot([baseline_times(1) baseline_times(1)],get(gca,'ylim'),'b--')
    plot([baseline_times(2) baseline_times(2)],get(gca,'ylim'),'b--')
    title(sprintf('%s %s ad times %1.1f %1.1f',file_name,ad_ch,ad_time,offsetTime))

    nexttile
    plot(chunk_times,chunks)
    hold on
    plot([ad_time ad_time],get(gca,'ylim'),'g--')
    plot([offsetTime offsetTime],get(gca,'ylim'),'r--')
    
    nexttile
    plot(on_chunk_times,on_chunks)
    hold on
    plot([ad_time ad_time],get(gca,'ylim'),'g--')
    plot([onsetTime onsetTime],get(gca,'ylim'),'b--')

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
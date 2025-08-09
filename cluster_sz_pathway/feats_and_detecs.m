%% plot_detections (with missed seizures first)

%% 1. Paths / env ---------------------------------------------
locations = seizure_termination_paths;
addpath(genpath(locations.annotation_grabber_folder))
addpath(genpath(locations.script_folder));
addpath(genpath(locations.ieeg_folder));
data_folder    = [locations.main_folder, 'data/'];
results_folder = [locations.main_folder, 'results/'];
out_folder     = [results_folder, 'cluster/'];

pwfile     = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

notchQ          = 10;
f0              = 60;
winSec = 1;   
lowBand  = [0 20];
highBand = [30 100];   
% -------------------------------------------------------------

% Info about bipolar channels
cluster_sz_data_file = '../../data/cluster_sz_data.xlsx';  
T = readtable(cluster_sz_data_file);

% ==== INPUT: CSVs ====
allDetectionsCSV    = '../../results/cluster/HUP199_phaseII_detections.csv';
falsePositivesCSV   = '../../results/cluster/HUP199_phaseII_false_positives.csv';
missedSeizuresCSV   = '../../results/cluster/HUP199_phaseII_missed_seizures.csv';

allDetectionsTable  = readtable(allDetectionsCSV);
falsePositivesTable = readtable(falsePositivesCSV);

% Try to read "missed" file if present
missedTimes = [];
if exist(missedSeizuresCSV, 'file')
    try
        missedTable = readtable(missedSeizuresCSV);
        missedTimes = get_time_column(missedTable);
    catch ME
        warning('⚠️ Could not read missed seizures table');
    end
else
    warning('⚠️ Missed seizures file not found: %s', missedSeizuresCSV);
end

% Assumes the detections table has a time-like column
allDetections   = get_time_column(allDetectionsTable);

% False detection times
falseDetections = get_time_column(falsePositivesTable);

% Extract filename from CSV name
[~, baseName, ~] = fileparts(allDetectionsCSV);  % e.g., 'HUP199_phaseII_detections'
fileName = erase(baseName, '_detections');       % 'HUP199_phaseII'

% Lookup channel info
rows = T(strcmp(T.filename, fileName), :);
if isempty(rows)
    error('No channel info found in cluster_sz_data.xlsx for %s', fileName);
end
ch1 = rows.bipolar_ch1{1};
ch2 = rows.bipolar_ch2{1};

%% Common params
preTime  = 15; % seconds before
postTime = 5;  % seconds after
tolerance = 0.01;  % seconds, for floating point match

figure; set(gcf,'position',[100 100 1300 650])

%% 2A. Loop over MISSED SEIZURES first ------------------------
for i = 1:numel(missedTimes)
    thisTime = missedTimes(i);
    t0 = thisTime - preTime;
    t1 = thisTime + postTime;

    fprintf('⏳ Processing MISSED seizure %d/%d at time %.3f s\n', i, numel(missedTimes), thisTime);

    % Download data
    try
        data = download_ieeg_data_sz(fileName, login_name, pwfile, [t0, t1], 1);
    catch ME
        warning('⚠️ Could not download data for %s at time %.3f s: %s', fileName, thisTime, ME.message);
        continue;
    end

    % Plot with "missed" styling
    plot_one_window(data, ch1, ch2, preTime, winSec, lowBand, highBand, ...
        sprintf('%s: MISSED seizure %d (%.3f s)', fileName, i, thisTime), ...
        'Missed Seizure', [0 0.4470 0.7410]); % MATLAB default blue

    disp('Press any key to continue...');
    pause
end

%% 2B. Loop over all DETECTIONS --------------------------------
for i = 1:length(allDetections)
    thisTime = allDetections(i);
    t0 = thisTime - preTime;
    t1 = thisTime + postTime;

    fprintf('⏳ Processing detection %d/%d at time %.3f s\n', i, length(allDetections), thisTime);

    % Download data
    try
        data = download_ieeg_data_sz(fileName, login_name, pwfile, [t0, t1], 1);
    catch ME
        warning('⚠️ Could not download data for %s at time %.3f s: %s', fileName, thisTime, ME.message);
        continue;
    end

    % Is this a false detection?
    isFalse  = any(abs(falseDetections - thisTime) < tolerance);
    labelStr = ternary(isFalse, 'False Detection', 'True Detection');
    lineCol  = ternary(isFalse, [1 0 0], [0 0.6 0]);  % red vs green

    % Plot with detection styling
    plot_one_window(data, ch1, ch2, preTime, winSec, lowBand, highBand, ...
        sprintf('%s: Detection %d (%.3f s)', fileName, i, thisTime), ...
        labelStr, lineCol);

    disp('Press any key to continue...');
    pause
end

%% ================= Helper functions =========================
function plot_one_window(data, ch1, ch2, preTime, winSec, lowBand, highBand, ttl, markerLabel, markerColor)
    fs        = data.fs;
    chLabels  = data.chLabels(:,1);
    raw       = data.values;

    % Notch filter design
    f0   = 60; notchQ = 10;
    wo = f0 / (fs/2);
    bw = wo / notchQ;
    [b, a] = iirnotch(wo, bw);

    % Find channel indices
    chanIdx1 = find(strcmp(chLabels, ch1));
    chanIdx2 = find(strcmp(chLabels, ch2));
    if isempty(chanIdx1) || isempty(chanIdx2)
        warning('⚠️ Channel(s) %s or %s not found in downloaded data.', ch1, ch2);
        return;
    end

    % Bipolar reference and notch
    bipolarData = raw(:, chanIdx2) - raw(:, chanIdx1);
    bipolarData = filtfilt(b, a, bipolarData);

    % Time vectors
    N   = size(bipolarData,1);
    t   = (0:N-1)/fs;
    detX = preTime; % event at center

    % --------- Compute features in 100 ms windows ----------
    winSamp   = round(winSec * fs);
    stepSamp  = winSamp;  % non-overlapping; change to smaller for overlap
    starts    = 1:stepSamp:(N - winSamp + 1);
    nWins     = numel(starts);

    LL        = zeros(nWins,1);
    ratio     = zeros(nWins,1);
    tWin      = zeros(nWins,1);

    for w = 1:nWins
        idx = starts(w):(starts(w)+winSamp-1);
        seg = bipolarData(idx);

        % Line length
        LL(w) = sum(abs(diff(seg)));

        % Power ratio (keep your current choice: using high-band power only)
        pLow  = bandpower(seg, fs, lowBand);
        pHigh = bandpower(seg, fs, [highBand(1) highBand(2)]);
        % ratio(w) = pHigh / (pLow + eps);   % original ratio if you want it
        ratio(w) = pHigh;                     % your current code

        % Time at the center of the window
        tWin(w) = (idx(1) + idx(end))/2 / fs;
    end

    % --------------------- Plotting ------------------------
    clf;
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    % 1) EEG
    ax1 = nexttile;
    plot(t, bipolarData, 'k');
    hold on;
    xline(ax1, detX, '--', markerLabel, 'Color', markerColor, ...
          'LabelVerticalAlignment','bottom','LineWidth',1.2);
    hold off;
    title(ttl);
    ylabel('Amplitude (\muV)');
    grid on;

    % 2) Line length
    ax2 = nexttile;
    plot(tWin, LL, 'LineWidth', 1);
    hold on;
    xline(ax2, detX, '--', 'Color', markerColor, 'LineWidth', 1.2);
    hold off;
    ylabel(sprintf('Line length (%.0f ms win)', winSec*1000));
    grid on;

    % 3) Power (your "ratio" series)
    ax3 = nexttile;
    plot(tWin, ratio, 'LineWidth', 1);
    hold on;
    xline(ax3, detX, '--', 'Color', markerColor, 'LineWidth', 1.2);
    hold off;
    xlabel('Time (s)');
    ylabel('High-band power');
    grid on;

    linkaxes([ax1 ax2 ax3], 'x');
    xlim([0 t(end)]);
end

function out = ternary(condition, trueVal, falseVal)
    if condition
        out = trueVal;
    else
        out = falseVal;
    end
end

function times = get_time_column(T)
    % Flexible extractor for a "time (sec)" column from a table.
    % Priority:
    % 1) Any varname containing 'time' (case-insensitive)
    % 2) Any varname ending with '_sec' or containing 'sec'
    % 3) First numeric column
    vn = string(T.Properties.VariableNames);

    % 1) contains 'time'
    idx = find(contains(lower(vn), 'time'), 1);
    if ~isempty(idx) && isnumeric(T.(vn(idx)))
        times = double(T.(vn(idx)));
        return;
    end

    % 2) endswith '_sec' or contains 'sec'
    idx = find(endsWith(lower(vn), '_sec') | contains(lower(vn), 'sec'), 1);
    if ~isempty(idx) && isnumeric(T.(vn(idx)))
        times = double(T.(vn(idx)));
        return;
    end

    % 3) first numeric column
    for k = 1:numel(vn)
        if isnumeric(T.(vn(k)))
            times = double(T.(vn(k)));
            return;
        end
    end

    error('No suitable time column found in table.');
end

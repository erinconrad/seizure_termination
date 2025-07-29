%% plot_detections

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

% Info about bipolar channels
cluster_sz_data_file = '../../data/cluster_sz_data.xlsx';  
T = readtable(cluster_sz_data_file);

% ==== INPUT: CSVs ====
allDetectionsCSV    = '../../results/cluster/HUP199_phaseII_detections.csv';
falsePositivesCSV   = '../../results/cluster/HUP199_phaseII_false_positives.csv';

allDetectionsTable  = readtable(allDetectionsCSV);
falsePositivesTable = readtable(falsePositivesCSV);

% Assumes both files have a column named "time" (in seconds)
allDetections  = allDetectionsTable.SeizureTime_sec;
falseDetections = falsePositivesTable.FalseDetection_sec;

% Extract filename from CSV name
[~, baseName, ~] = fileparts(allDetectionsCSV);  % e.g., 'HUP192_phaseII_detections'
fileName = erase(baseName, '_detections');       % 'HUP192_phaseII'

% Lookup channel info
rows = T(strcmp(T.filename, fileName), :);
if isempty(rows)
    error('No channel info found in cluster_sz_data.xlsx for %s', fileName);
end
ch1 = rows.bipolar_ch1{1};
ch2 = rows.bipolar_ch2{1};

%% 2. Loop over all detections ---------------------------------------------
preTime = 15; % seconds before
postTime = 5; % seconds after
tolerance = 0.01;  % seconds, for floating point match

figure;
for i = 10:length(allDetections)
    thisTime = allDetections(i);
    t0 = thisTime - preTime;
    t1 = thisTime + postTime;

    fprintf('⏳ Processing detection %d/%d at time %.1f s\n', i, length(allDetections), thisTime);

    % Download 20s of data around the detection
    try
        data = download_ieeg_data_sz(fileName, login_name, pwfile, [t0, t1], 1);
    catch ME
        warning('⚠️ Could not download data for %s at time %.1f s: %s', fileName, thisTime, ME.message);
        continue;
    end

    fs        = data.fs;
    chLabels  = data.chLabels(:,1);
    raw       = data.values;

    % Design notch filter
    wo = f0 / (fs/2);     % Normalize frequency
    bw = wo / notchQ;          % Bandwidth
    [b, a] = iirnotch(wo, bw);  % Design notch filter

    % Find channel indices
    chanIdx1 = find(strcmp(chLabels, ch1));
    chanIdx2 = find(strcmp(chLabels, ch2));
    if isempty(chanIdx1) || isempty(chanIdx2)
        warning('⚠️ Channel(s) %s or %s not found in downloaded data for %s.', ch1, ch2, fileName);
        continue;
    end

    % Bipolar reference
    bipolarData = raw(:, chanIdx2) - raw(:, chanIdx1);

    % apply a notch filter
    bipolarData = filtfilt(b, a, bipolarData );

    % Determine if this is a false detection
    isFalse = any(abs(falseDetections - thisTime) < tolerance);
    labelStr = ternary(isFalse, 'False Detection', 'True Detection');
    lineColor = ternary(isFalse, 'r--', 'g--');

    % Plot
    t = linspace(0, 20, size(bipolarData, 1));
    
    plot(t, bipolarData);
    hold on;
    xline(15, lineColor, labelStr, 'LabelVerticalAlignment', 'bottom');
    hold off;

    title(sprintf('%s: Detection %d (%.1f s)', fileName, i, thisTime));
    xlabel('Time (s)');
    ylabel('Amplitude (uV)');
    grid on;

    disp('Press any key to continue to the next detection...');
    pause
    
end

%% Helper function (ternary)
function out = ternary(condition, trueVal, falseVal)
    if condition
        out = trueVal;
    else
        out = falseVal;
    end
end

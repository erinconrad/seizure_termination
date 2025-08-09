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
winSec = 0.50;   % 100 ms windows
lowBand  = [0 20];
highBand = [20 50];   
% -------------------------------------------------------------

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
set(gcf,'position',[100 100 1300 650])
for i = 1:length(allDetections)
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
    wo = f0 / (fs/2);    % Normalize frequency
    bw = wo / notchQ;    % Bandwidth
    [b, a] = iirnotch(wo, bw);

    % Find channel indices
    chanIdx1 = find(strcmp(chLabels, ch1));
    chanIdx2 = find(strcmp(chLabels, ch2));
    if isempty(chanIdx1) || isempty(chanIdx2)
        warning('⚠️ Channel(s) %s or %s not found in downloaded data for %s.', ch1, ch2, fileName);
        continue;
    end

    % Bipolar reference and notch
    bipolarData = raw(:, chanIdx2) - raw(:, chanIdx1);
    bipolarData = filtfilt(b, a, bipolarData);

    % Time vectors
    N   = size(bipolarData,1);
    t   = (0:N-1)/fs;                 % 0–20 s
    detX = preTime;                   % detection marker at 15 s into the 20 s window

    % Is this a false detection?
    isFalse  = any(abs(falseDetections - thisTime) < tolerance);
    labelStr = ternary(isFalse, 'False Detection', 'True Detection');
    lineCol  = ternary(isFalse, [1 0 0], [0 0.6 0]);  % red vs green

    %% --------- Compute features in 100 ms windows ----------
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

        % Power ratio
        pLow  = bandpower(seg, fs, lowBand);
        pHigh = bandpower(seg, fs, [highBand(1) highBand(2)]);
        %ratio(w) = pHigh / (pLow + eps);
        ratio(w) = pHigh;

        % Time at the center of the window
        tWin(w) = (idx(1) + idx(end))/2 / fs;
    end

    %% --------------------- Plotting ------------------------
    clf;  % clear the figure each iteration
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    % 1) EEG
    ax1 = nexttile;
    plot(t, bipolarData, 'k');
    hold on;
    xline(ax1, detX, '--', labelStr, 'Color', lineCol, ...
          'LabelVerticalAlignment','bottom','LineWidth',1.2);
    hold off;
    title(sprintf('%s: Detection %d (%.1f s)', fileName, i, thisTime));
    ylabel('Amplitude (\muV)');
    grid on;

    % 2) Line length
    ax2 = nexttile;
    plot(tWin, LL, 'LineWidth', 1);
    hold on;
    xline(ax2, detX, '--', 'Color', lineCol, 'LineWidth', 1.2);
    hold off;
    ylabel(sprintf('Line length (%.0f ms win)', winSec*1000));
    grid on;

    % 3) Power ratio >15 Hz / <15 Hz
    ax3 = nexttile;
    plot(tWin, ratio, 'LineWidth', 1);
    hold on;
    xline(ax3, detX, '--', 'Color', lineCol, 'LineWidth', 1.2);
    hold off;
    xlabel('Time (s)');
    ylabel('Power ratio >15 / <15 Hz');
    grid on;

    linkaxes([ax1 ax2 ax3], 'x');
    xlim([0 t(end)]);

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
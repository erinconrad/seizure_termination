%% measure_accuracy_all.m
%% CONFIG -----------------------------------------------------------------
seizureSpreadsheet = '../../results/cluster/all_seizure_times.csv';   % includes all true seizure times
detectFolder       = '../../results/cluster/';                        % where detection CSVs live
matchWindow        = 120;   % seconds (±) for a valid detection
timeColTrue        = 'start';   % column in spreadsheet with seizure time
filenameCol        = 'IEEGname'; % column with filenames (match detection CSVs)

% Define per‑file search intervals in seconds
searchIntervals = struct( ...
    'HUP230_phaseII', [3860.24-400,257568.67+100],...
    'HUP247_phaseII', [60052.43-400, 80786.43+100], ...
    'HUP251_phaseII', [30043.37207-400, 35440.62988+100], ...
    'HUP266_phaseII', [253184.72-400, 301581.35+100],...
    'HUP241_phaseII', [41248.98-400, 171393.73+100],...
    'HUP213_phaseII_D01',[14988.18-400,65008.3301],...
    'HUP259_phaseII',[586042.29-400,776658.81+100]...
);
% -------------------------------------------------------------------------

%% Load ground truth seizure times
trueTbl = readtable(seizureSpreadsheet);

% Get unique filenames
uniqueFiles = fieldnames(searchIntervals);

% Summary accumulators
all_stats = [];

for f = 1:numel(uniqueFiles)
    fname = uniqueFiles{f};
    fprintf('\n=== Evaluating %s ===\n', fname);

    % Get true seizure times for this file and restrict to search interval
    interval = searchIntervals.(fname);  % [start, end] in seconds
    isThisFile = strcmp(trueTbl.(filenameCol), fname);
    allSeiz    = double(trueTbl.(timeColTrue)(isThisFile));
    seizTimes  = sort(allSeiz(allSeiz >= interval(1) & allSeiz <= interval(2)));

    % Load detections
    detFile = fullfile(detectFolder, [fname '_detections.csv']);
    if ~isfile(detFile)
        fprintf('  ⚠ Detection file not found: %s\n', detFile);
        continue
    end
    detTbl = readtable(detFile);
    detTimes = sort(double(detTbl.SeizureTime_sec));

    % Initialize match tracking
    latencies = [];
    usedDetIdx = false(size(detTimes));

    % Match each seizure to closest detection (within window)
    for sz = seizTimes(:).'
        deltas = abs(detTimes - sz);
        [dmin, bestIdx] = min(deltas);
        if dmin <= matchWindow
            latencies(end+1,1) = detTimes(bestIdx) - sz; %#ok<AGROW>
            usedDetIdx(bestIdx) = true;
        end
    end

    % Metrics
    nSeiz      = numel(seizTimes);
    nDetected  = numel(latencies);
    nFalse     = sum(~usedDetIdx);

    interval = searchIntervals.(fname);
    duration = (interval(2) - interval(1)) / 3600; % in hours

    sensitivity      = nDetected / nSeiz;
    falseRatePerHour = nFalse / duration;

    % Print summary
    fprintf('  Seizures:                  %d\n', nSeiz);
    fprintf('  Detected:                  %d\n', nDetected);
    fprintf('  Sensitivity:               %.2f %%\n', 100 * sensitivity);
    fprintf('  False detections:          %d\n', nFalse);
    fprintf('  Recording duration (hrs):  %.2f\n', duration);
    fprintf('  False detections / hr:     %.2f\n', falseRatePerHour);

    % Optionally plot latency histogram
    figure;
    edges = -matchWindow:1:matchWindow;
    histogram(latencies, edges, ...
        'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
    xline(0, '--', 'LineWidth', 1);
    xlabel('Detection time – Seizure time (s)');
    ylabel('Count');
    title(sprintf('Latency distribution (%s)', fname));
    grid on;

    % Save per-file stats
    all_stats = [all_stats; table({fname}, nSeiz, nDetected, sensitivity, ...
        nFalse, duration, falseRatePerHour, ...
        'VariableNames', {'Filename', 'nSeizures', 'nDetected', 'Sensitivity', ...
                          'nFalsePositives', 'DurationHours', 'FalsePosRatePerHour'})];

    %% Save false positives
    falsePosTimes = detTimes(~usedDetIdx);
    falsePosTable = table(falsePosTimes, 'VariableNames', {'FalseDetection_sec'});
    fp_outfile = fullfile(detectFolder, [fname '_false_positives.csv']);
    writetable(falsePosTable, fp_outfile);
    fprintf('  ➜ Saved %d false positives to %s\n', numel(falsePosTimes), fp_outfile);

    %% Save missed seizures
    matchedSeizureIdx = false(size(seizTimes));  % same size as seizTimes
    for i = 1:length(seizTimes)
        sz = seizTimes(i);
        deltas = abs(detTimes - sz);
        if any(deltas <= matchWindow)
            matchedSeizureIdx(i) = true;
        end
    end
    missedSeizures = seizTimes(~matchedSeizureIdx);
    missedTable = table(missedSeizures, 'VariableNames', {'MissedSeizure_sec'});
    missed_outfile = fullfile(detectFolder, [fname '_missed_seizures.csv']);
    writetable(missedTable, missed_outfile);
    fprintf('  ➜ Saved %d missed seizures to %s\n', numel(missedSeizures), missed_outfile);

end

% Save summary table
summaryOut = fullfile(detectFolder, 'detection_accuracy_summary.csv');
writetable(all_stats, summaryOut);
fprintf('\nSaved summary to %s\n', summaryOut);

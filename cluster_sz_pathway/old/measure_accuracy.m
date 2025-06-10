%% measure_accuracy
%% CONFIG -----------------------------------------------------------------
seizureCsv   = '../../results/cluster/HUP247_phaseII_true.csv';      % reference seizure times
detectCsv    = '../../results/cluster/HUP247_phaseII_detections.csv';% detector output
timeVar      = 'SeizureTime_sec';    % column holding times (s from start)
matchWindow  = 20;                   % hit window (± s)
% -------------------------------------------------------------------------

%% 1) Load the two CSVs
seizTbl = readtable(seizureCsv);
detTbl  = readtable(detectCsv);

% numeric, column vectors, sorted
seizTimes = sort(double(seizTbl.(timeVar)));
detTimes  = sort(double(detTbl.(timeVar)));

%% 2) Match each seizure to the nearest detection within matchWindow
latencies  = [];                           % detection – seizure (signed)
usedDetIdx = false(size(detTimes));        % detections already claimed

for sz = seizTimes(:).'                    % loop seizures
    deltas = abs(detTimes - sz);
    [dmin, bestIdx] = min(deltas);
    if dmin <= matchWindow
        latencies(end+1,1) = detTimes(bestIdx) - sz; %#ok<AGROW>
        usedDetIdx(bestIdx) = true;         % prevent re‑use
    end
end

%% 3) Summary statistics
nSeizures         = numel(seizTimes);
nDetectedSeizures = numel(latencies);
sensitivity       = nDetectedSeizures / nSeizures;

nFalseDetections  = sum(~usedDetIdx);

% duration (h) – swap in explicit start/stop if you have them
recordingHours    = (max([seizTimes; detTimes]) - min([seizTimes; detTimes])) / 3600;
falseRatePerHour  = nFalseDetections / recordingHours;

fprintf('Sensitivity                 : %.2f %%\n', 100*sensitivity);
fprintf('False detections            : %d\n', nFalseDetections);
fprintf('Recording duration (hours)  : %.2f\n', recordingHours);
fprintf('False detections per hour   : %.2f\n', falseRatePerHour);

%% 4) Histogram of latencies for true positives
figure;
edges = -matchWindow:1:matchWindow;      % 1‑second bins
histogram(latencies, edges, ...
    'FaceColor',[0.25 0.45 0.85], 'EdgeColor','none');
xline(0,'--','LineWidth',1);
xlabel('Detection time – Seizure time (s)');
ylabel('Count');
title('Latency distribution for true detections');
grid on;

%% 5) Print the false detections  %%% NEW ---------------------------------
falseDetTimes = detTimes(~usedDetIdx);
if isempty(falseDetTimes)
    fprintf('No false detections were found.\n');
else
    fprintf('\nFalse detection times (s):\n');
    fprintf('  %.2f\n', falseDetTimes);   % each on its own line
end

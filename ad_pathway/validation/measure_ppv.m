%% measure_ppv

clear;
clc;

%% Paths
locations = seizure_termination_paths;
addpath(genpath(locations.annotation_grabber_folder))
data_folder = [locations.main_folder,'data/'];
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'detections/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Get list of detection files
listing = dir([out_folder,'*detections.csv']);
nfiles = length(listing);

n_ads = zeros(nfiles,1);
ppv = nan(nfiles,1);

% Loop over files
for i = 1:nfiles
    if strcmp(listing(i).name(1),'.'), continue; end

    fname = [out_folder,'/',listing(i).name];
    currT = readtable(fname);

    % Loop over ads
    ad_rows = find(strcmp(currT.Type,'AD'));
    if isempty(ad_rows), continue; end

    % get my designations
    erin_des = currT.Erin(ad_rows);
    curr_ppv = sum(strcmp(erin_des,'y'))/length(ad_rows);
    ppv(i) = curr_ppv;
    n_ads(i) = length(ad_rows);

end

%% Summarize
fprintf(['\nThe median (IQR) number of ad detections was %1.1f (%1.1f-%1.1f).\n'],...
    median(n_ads),prctile(n_ads,25),prctile(n_ads,75));

fprintf(['\nThe median (IQR) PPV of ad detections was %1.1f (%1.1f-%1.1f).\n'],...
    nanmedian(ppv),prctile(ppv,25),prctile(ppv,75));

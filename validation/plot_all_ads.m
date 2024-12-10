%% plot_all_ad
%% Parameters
surr_time = 15;

%% Paths
locations = seizure_termination_paths;
addpath(genpath(locations.annotation_grabber_folder))
data_folder = [locations.main_folder,'data/'];
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'detections/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% ieeg stuff
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

mkdir(out_folder)


%% Get list of detection files
listing = dir([out_folder,'*.detections.csv']);
nfiles = length(listing);

% Loop
for i = 1:nfiles
    fname = [out_folder,'/',listing(i).name];
    currT = readtable(fname);

    %% PLot
    plot_ad_detections(currT,out_folder,surr_time)

end
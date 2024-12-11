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
listing = dir([out_folder,'*detections.csv']);
nfiles = length(listing);

% Loop
for i = 1:nfiles
    if strcmp(listing(i).name(1),'.'), continue; end

    fname = [out_folder,'/',listing(i).name];
    currT = readtable(fname);
    ieeg_name = currT.FileName{1};
    modifier = currT.Modifier{1};
    ad_out_name = [out_folder,ieeg_name,'_',sprintf('%d',modifier),'_ads'];


    %% PLot
    plot_ad_detections(currT,ad_out_name,surr_time)

end
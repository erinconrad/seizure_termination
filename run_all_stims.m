% do all stim

%% Parameters
overwrite = 1;
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

%% Get list of stim sessions
fT = readtable([results_folder,'hfs_sessions.csv']);
nfiles = size(fT,1);

%% Initialize variables
% initialize tracker of ad and stim
T = table('Size',[0 5],'VariableTypes',{'cell','cell','cell','double','double'},'VariableNames',{'FileName','Modifier','Type','Channels','OnTime','OffTime'});


% Loop over files
for i = 1:nfiles
    ieeg_name = fT.ieeg_name{i};
    start_time = fT.start(i);
    end_time = fT.xEnd(i);
    modifier = fT.Modifier(i);

    if overwrite == 0
        if exist([out_folder,ieeg_name,'_',modifier,'_detections.csv'],'file') ~= 0
            fprintf('\nAlready did %s, skipping...\n',ieeg_name);
            continue;
        end
    end

    %% Get the stim and AD annotations
    fprintf('\nDoing %s_%d, file %d of %d\n',ieeg_name,modifier,i,nfiles);
    currT = find_ad_fcn(ieeg_name,start_time,end_time);
    currT.FileName = repmat({ieeg_name},size(currT,1),1);
    currT.Modifier = repmat(modifier,size(currT,1),1);
    currT = currT(:, ['FileName', 'Modifier', T.Properties.VariableNames(1:end-2)]);
    
    writetable(currT,[out_folder,ieeg_name,'_',modifier,'_detections.csv']);

    %% Plot the EEG surrounding each of the ADs
    plot_ad_detections(currT,out_folder,surr_time)
end


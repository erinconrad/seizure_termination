% do all stim

%% Parameters
overwrite = 1;
surr_time = 15;
max_dur = 2000;

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



% Loop over files
for i = 1:nfiles
    ieeg_name = fT.ieeg_name{i};
    start_time = fT.start(i);
    end_time = fT.xEnd(i);
    modifier = fT.Modifier(i);
    csv_out_name = [out_folder,ieeg_name,'_',sprintf('%d',modifier),'_detections.csv'];

    if overwrite == 0
        if exist(csv_out_name,'file') ~= 0
            fprintf('\nAlready did %s modifier %d, skipping...\n',ieeg_name,modifier);
            continue;
        end
    end

    %% Get the stim and AD annotations
    fprintf('\nDoing %s_%d, file %d of %d...',ieeg_name,modifier,i,nfiles);

    %% Initialize variables
    % initialize tracker of ad and stim
    T = table('Size',[0 7],'VariableTypes',{'cell','double','double','cell','cell','double','double'},...
        'VariableNames',{'FileName','Modifier','Run','Type','Channels','OnTime','OffTime'});


    %% Handle long stim sessions
    duration = end_time - start_time;
    nruns = ceil(duration/max_dur);
    for r = 1:nruns
        fprintf('    Doing run %d of %d\n',r,nruns);
        curr_start = start_time + max_dur*(r-1);
        curr_end = min([end_time,curr_start+max_dur]);

        tempT = find_ad_fcn(ieeg_name,curr_start,curr_end);
        currT = tempT;
        currT.FileName = repmat({ieeg_name},size(currT,1),1);
        currT.Modifier = repmat(modifier,size(currT,1),1);
        currT.Run = repmat(r,size(currT,1),1);
        currT = currT(:, ['FileName', 'Modifier','Run', currT.Properties.VariableNames(1:end-3)]);
        T = [T;currT];
    end

    
    
    writetable(T,csv_out_name);

    %% Plot the EEG surrounding each of the ADs
    plot_ad_detections(T,out_folder,surr_time)
end


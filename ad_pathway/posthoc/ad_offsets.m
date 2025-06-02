%% add_offsets

%{
This script loops over a directory of detections, and adds the offset to
each detection

%}

%% Paths
locations = seizure_termination_paths;
addpath(genpath(locations.annotation_grabber_folder))
data_folder = [locations.main_folder,'data/'];
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'detections_test/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% ieeg stuff
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

% Loop over the detection files
listing = dir([out_folder,'*detections.csv']);
for i = 1:length(listing)
    
    curr_file = [out_folder,listing(i).name];

    % load it
    currT = readtable(curr_file);

    % get offset times for each AD detection
    ad_rows = find(strcmp(currT.Type,'AD'));
    for j = 1:length(ad_rows)
        row = ad_rows(j);
        file_name = currT.FileName{row};
        ad_time = currT.OnTime(row);
        
        % get the most recent stim time
        earlier_stim_rows = find(strcmp(currT.Type(1:row),'stim'));
        last_stim_row = earlier_stim_rows(end);
        if isempty(last_stim_row)
            error('why no stim??')
        end
        stim_times = [currT.OnTime(last_stim_row) currT.OffTime(last_stim_row)];
        ad_ch = currT.Channels{row};

        % get the offset time
        offset = get_precise_offset(file_name,ad_time,ad_ch,stim_times);

        currT.OffTime(row) = offset;
        
    end

    % rewrite the file with this offset information
    writetable(currT,curr_file);

    ad_out_name = strrep(curr_file,'detections.csv','ads');

    % Plot the EEG surrounding each of the ADs
    surr_time = 15;
    plot_ad_detections(currT,ad_out_name,surr_time)
end
% do all stim

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

%% Get list of stim sessions
fT = readtable([results_folder,'hfs_sessions_fake.csv']);
nfiles = size(fT,1);

%% Initialize variables
% initialize tracker of ad and stim
T = table('Size',[0 5],'VariableTypes',{'cell','cell','cell','double','double'},'VariableNames',{'FileName','Type','Channels','OnTime','OffTime'});


% Loop over files
for i = 1:nfiles
    ieeg_name = fT.ieeg_name{i};
    start_time = fT.start(i);
    end_time = fT.xEnd(i);

    %% Get the stim and AD annotations
    fprintf('\nDoing %s, file %d of %d\n',ieeg_name,i,nfiles);
    currT = find_ad_fcn(ieeg_name,start_time,end_time);
    currT.FileName = repmat({ieeg_name},size(currT,1),1);
    currT = currT(:, ['FileName', T.Properties.VariableNames(1:end-1)]);
    
    writetable(currT,[out_folder,ieeg_name,'_detections.csv']);

    %% Plot the EEG surrounding each of the ADs
    ad_rows = find(strcmp(currT.Type,'AD'));
    if isempty(ad_rows), continue; end
    figure
    set(gcf,'position',[1 1 1400 200*length(ad_rows)])
    tiledlayout(length(ad_rows),1,'tilespacing','tight','padding','tight')
    for ia = 1:length(ad_rows)
        nexttile
        row = ad_rows(ia);
        ad_time = currT.OnTime(row);
        ad_ch = currT.Channels{row};
        data = download_ieeg_data(ieeg_name, login_name, pwfile, ...
            [ad_time-surr_time,ad_time+surr_time], 1);
        chLabels = data.chLabels(:,1);
        chLabels = decompose_labels(chLabels);
        ch_idx = strcmp(chLabels,ad_ch);
        values = data.values(:,ch_idx);
        plot(linspace(-surr_time,surr_time,length(values)),values)
        xlabel('seconds')
        title(sprintf('%1.1f %s',ad_time,ad_ch))
        
    end
    print(gcf,[out_folder,ieeg_name,'_ads'],'-dpng')
    close gcf

end


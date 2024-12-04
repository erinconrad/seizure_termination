%% find HFS stim sessions

%% Parameters
rows = 16:59;

%% Paths
locations = seizure_termination_paths;
addpath(genpath(locations.annotation_grabber_folder))
data_folder = [locations.main_folder,'data/'];
results_folder = [locations.main_folder,'results/'];

%% Grab file names of interest
dT = readtable([data_folder,'master_pt_list.xlsx']);
hup_nos = dT.HUPID;
poss_fnames = cellfun(@(x) [x,'_phaseII'], hup_nos,'uniformoutput',false);
filenames = poss_fnames(rows);
nfiles = length(filenames);

%% Get all annotations
% initialize empty table
aT = table('Size',[0,3],...
    'VariableTypes',{'string','double','string'},...
    'VariableNames',{'ieeg_name','annotation_times','annotations'});

for i = 1:nfiles
    ieeg_name = filenames{i};

    %% Get the annotations
    try
        currT = grab_all_annotations(ieeg_name);
    catch
        fprintf('\nfailed %s\n',ieeg_name);
        continue;
    end

    if isempty(currT)
        continue;
    end

    %% Mine for stim-y annotations
    annotations = currT.annotations;
    stim_annotations = contains(annotations,["stim","hz"],'IgnoreCase',true);
    curr_stimT = currT(stim_annotations,:);

    curr_stimT.ieeg_name = repmat({ieeg_name},size(curr_stimT,1),1);
    curr_stimT = curr_stimT(:, ['ieeg_name', curr_stimT.Properties.VariableNames(1:end-1)]);
    aT = [aT;curr_stimT];

    %% Save
    writetable(aT,[results_folder,'poss_stim.csv'])
end




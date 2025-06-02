%% gui_validation

clear;
clc;

%% Parameters
surr_time = 9;

% overwrite = 0 (default): skip previously completed files
% overwrite = 1: re-write selections of all files
% overwrite = 2: show all files, but do not save any selections
overwrite = 0;
gain = 0;

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

%{ ----- NO NEED TO EDIT BELOW THIS ------------ %}

%% Get list of detection files
listing = dir([out_folder,'*detections.csv']);
nfiles = length(listing);

% Loop over files
for i = 1:nfiles
    if strcmp(listing(i).name(1),'.'), continue; end

    fname = [out_folder,'/',listing(i).name];
    currT = readtable(fname);
    ieeg_name = currT.FileName{1};
    modifier = currT.Modifier(1);

    fprintf('\nDoing %s\n',ieeg_name);


    % See if it has a variable name for my designations, if not add it
    if sum(strcmp(currT.Properties.VariableNames,'Erin')) == 0
        currT.Erin = repmat({''},size(currT,1),1);
    end

    % Loop over ads
    ad_rows = find(strcmp(currT.Type,'AD'));
    if isempty(ad_rows), continue; end

    for j = 1:length(ad_rows)
        r = ad_rows(j);
        if overwrite == 0 
            if strcmp(currT.Erin{r},'y') || strcmp(currT.Erin{r},'n')
                continue
            end
        end

        % Get AD details
        ad_time = currT.OnTime(r);
        ad_ch = currT.Channels{r};
        ieeg_name = currT.FileName{r};

        % Download the data
        data = download_ieeg_data(ieeg_name, login_name, pwfile, ...
            [ad_time-surr_time,ad_time+surr_time], 1);
        chLabels = data.chLabels(:,1);
        chLabels = decompose_labels(chLabels);
        [~,bipolarIndices,altBipolarIndices] = find_bipolar_pairs(chLabels,1:length(chLabels));
        ch_idx = strcmp(chLabels,ad_ch);
        values = data.values;
        fs = data.fs;

        % bipolar
        biValues = values;
        biValues(:,~isnan(altBipolarIndices)) = values(:,~isnan(altBipolarIndices)) - values(:,altBipolarIndices(~isnan(altBipolarIndices)));
        biValues(:,isnan(altBipolarIndices)) = nan;
        biValues = biValues - mean(biValues,1,"omitnan");

        % restrict to ch of interest
        biValues = biValues(:,ch_idx);
    
        % notch filter to aid with visualization
        biValues = bandstop(biValues,[58 62],data.fs);
        
        % call this plot values
        plot_values = biValues;
        plot_labels = chLabels(ch_idx);

        % Create a figure for displaying the EEG data
        figure;
        hFig = gcf;
        set(hFig, 'Position', [10 10 1400 700]);
        tiledlayout(1, 1, 'Padding', 'compact');
        nexttile;
        
        % Initialize user_input
        user_input = '';  % Initialize as empty

        % Plot EEG data with prompt and selection texts
        prompt_text = plot_eeg_with_text(plot_values, fs, plot_labels, gain);
    

         %% Store Handles and Data in appdata
        setappdata(hFig, 'prompt_text', prompt_text);
        setappdata(hFig, 'gain', gain);
        setappdata(hFig, 'fs', fs);
        setappdata(hFig, 'plot_values', plot_values);
        setappdata(hFig, 'plot_labels', plot_labels);

        %% Handle User Input
        % Initialize variables and store them in appdata
        arrow_key_input = '';
        setappdata(hFig, 'user_input', user_input);
        setappdata(hFig, 'arrow_key_input', arrow_key_input);
       
        % Set up the figure to detect key presses using WindowKeyPressFcn
        set(hFig, 'WindowKeyPressFcn', @keyPressHandler);
        
        % Wait until the user has confirmed their response
        uiwait(hFig);  % Wait until uiresume is called

        %% Retrieve User Input and Other Data
        % Retrieve the stored arrow_key_input and user_input from appdata
        arrow_key_input = getappdata(hFig, 'arrow_key_input');
        user_input = getappdata(hFig, 'user_input');
        gain = getappdata(hFig, 'gain');
        
        %% Handle End of File and Write Responses
        % Now you can safely close the figure
        if ishandle(hFig)
            close(hFig);
        end

        % Store the filename and response in the responses array
        currT.Erin{r} = user_input;
        
        if overwrite ~= 2
            % Save the table as a CSV file (save as you go)
            writetable(currT, fname);
        end


    end


end



function prompt_text = plot_eeg_with_text(plot_values, fs, plot_labels, gain)
    % This function plots the EEG data and adds prompt and selection texts
    plot_scalp_eeg_clean(plot_values, fs, plot_labels, gain);
    
    % Add prompt text as super title
    prompt_str = 'Designate whether true AD (y) or not (n)';
    sgtitle(prompt_str, 'FontSize', 16, 'Color', 'blue');
    

    % Since sgtitle does not return a handle, set prompt_text to empty
    prompt_text = [];
end

function keyPressHandler(hObject, event)
    % Nested key press handler function

    % Retrieve data stored in appdata
    arrow_key_input = getappdata(hObject, 'arrow_key_input');
    gain = getappdata(hObject, 'gain');
    user_input = getappdata(hObject, 'user_input');
    plot_values = getappdata(hObject, 'plot_values');
    plot_labels = getappdata(hObject, 'plot_labels');
    fs = getappdata(hObject, 'fs');

    % Handle key presses
    if strcmp(event.Key, 'uparrow')
        disp('Up arrow pressed');
        arrow_key_input = 'up';
        gain = gain - 100;
        setappdata(hObject, 'gain', gain);
        hold off;
        
        % Retrieve current user_input for displaying
        current_input = getappdata(hObject, 'user_input');
        % Re-plot EEG with updated gain and existing user_input
        plot_eeg_with_text(plot_values, fs, plot_labels, gain);

        % Update appdata with new text handles
        setappdata(hObject, 'prompt_text', []); % sgtitle doesn't return a handle

    elseif strcmp(event.Key, 'downarrow')
        disp('Down arrow pressed');
        arrow_key_input = 'down';
        gain = gain + 100;
        setappdata(hObject, 'gain', gain);
        hold off;
        
        % Retrieve current user_input for displaying
        current_input = getappdata(hObject, 'user_input');
        % Re-plot EEG with updated gain and existing user_input
        plot_eeg_with_text(plot_values, fs, plot_labels, gain);

        % Update appdata with new text handles
        setappdata(hObject, 'prompt_text', []); % sgtitle doesn't return a handle
    elseif strcmp(event.Key, 'return')
        % User pressed Enter to confirm input
        if ismember(user_input, {'y', 'n'})
            disp(['User input confirmed: ', user_input]);
            % Store the user input
            setappdata(hObject, 'user_input', user_input);
            % Resume execution (but do NOT close the figure here)
            uiresume(hObject);
        else
            disp('Please press "y" (yes), or "n" (no) before pressing Enter to confirm.');
        end
    elseif ismember(event.Key, {'y', 'n'})
        % Store the user's input character but wait for Enter to confirm
        user_input = event.Key;
        disp(['You pressed "', user_input, '". Press Enter to confirm.']);
        % Store the tentative user input
        setappdata(hObject, 'user_input', user_input);
  
    else
        disp(['Other key pressed: ', event.Key]);
    end

   
end



function plot_scalp_eeg_clean(values, fs, labels, gain)
    % This function plots the EEG data in a scalp montage
   
    
    dur = size(values,1)/fs;
   
    
    
    plot(linspace(0, dur, size(values,1)), values, 'k');
    hold on;
    
    text(dur + 0.05, median(values), sprintf('%s', labels{1}), 'FontSize', 15);
    
    plot([dur/2 dur/2], get(gca, 'ylim'), 'r--');
    
    yl = get(gca,'ylim');
    yl = [yl(1)-gain yl(2) + gain];
    set(gca,'ylim',yl);

    yticklabels([]);
    xticks(1:floor(dur));
    xlabel('Time (seconds)');
    set(gca, 'FontSize', 15);
    
    %% Add second markers
    for i = 1:floor(dur)
        plot([i, i], get(gca, 'ylim'), 'k--');
    end
end

function [bipolarPairs,bipolarIndices,altBipolarIndices] = find_bipolar_pairs(contacts,contactIndices)

% Initialize cell array to store bipolar pairs
bipolarPairs = {};
bipolarIndices = [];
altBipolarIndices = nan(length(contacts),1);

% Loop through each contact to extract letter and number parts
for i = 1:numel(contacts)
    % Parse the current contact
    match = regexp(contacts{i}, '([A-Z]+)(\d+)', 'tokens');
    if ~isempty(match)
        letterPart = match{1}{1};               % Extract letter part
        numberPart = str2double(match{1}{2});    % Extract number part as double

        % Check for the next contact with the same letter part and incremented number
        nextContact = sprintf('%s%d', letterPart, numberPart + 1);
        if any(strcmp(contacts, nextContact))
            % If next contact exists, add the pair to bipolarPairs
            bipolarPairs = [bipolarPairs; {contacts{i}, nextContact}];
            bipolarIndices = [bipolarIndices; contactIndices(i) contactIndices(find(strcmp(contacts,nextContact)))];
            altBipolarIndices(i) = find(strcmp(contacts,nextContact));
        end
    end
end
end
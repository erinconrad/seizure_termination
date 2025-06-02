function plot_ad_detections(T,ad_out_name,surr_time)

num_rows = 5;

%% Paths
locations = seizure_termination_paths;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% ieeg stuff
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

currT = T;
ad_rows = find(strcmp(currT.Type,'AD'));
if isempty(ad_rows), return; end

nads = length(ad_rows);
nfigs = ceil(nads/num_rows);

for i = 1:nfigs

    start_ad = num_rows * (i-1) + 1;
    end_ad = min([num_rows*i,nads]);
    curr_rows = ad_rows(start_ad:end_ad);
    n_curr_ads = length(curr_rows);

    figure
    set(gcf,'position',[1 1 1400 200*length(curr_rows)])
    tiledlayout(n_curr_ads,1,'tilespacing','tight','padding','tight')
    for ia = 1:n_curr_ads
        nexttile
        row = curr_rows(ia);
        ad_time = currT.OnTime(row);
        ad_off_time = currT.OffTime(row);
        rel_off = ad_off_time - ad_time;
        ad_ch = currT.Channels{row};
        ieeg_name = currT.FileName{row};
        data = download_ieeg_data(ieeg_name, login_name, pwfile, ...
            [ad_time-surr_time,ad_time+surr_time], 1);
        chLabels = data.chLabels(:,1);
        chLabels = decompose_labels(chLabels);
        [~,bipolarIndices,altBipolarIndices] = find_bipolar_pairs(chLabels,1:length(chLabels));
        ch_idx = strcmp(chLabels,ad_ch);
        %values = data.values(:,ch_idx);
        values = data.values;

        % bipolar
        biValues = values;
        biValues(:,~isnan(altBipolarIndices)) = values(:,~isnan(altBipolarIndices)) - values(:,altBipolarIndices(~isnan(altBipolarIndices)));
        biValues(:,isnan(altBipolarIndices)) = nan;
        biValues = biValues - mean(biValues,1,"omitnan");

        % restrict to ch of interest
        biValues = biValues(:,ch_idx);
    
        % notch filter to aid with visualization
        biValues = bandstop(biValues,[58 62],data.fs);

        

        plot(linspace(-surr_time,surr_time,length(biValues)),biValues)

        hold on

        plot([rel_off rel_off],get(gca,'ylim'),'r--')
        xlim([-surr_time,surr_time])
        xlabel('seconds')
        title(sprintf('%1.1f %s',ad_time,ad_ch))
        
    end
    print(gcf,[ad_out_name,sprintf('_fig%d',i)],'-dpng')
    close gcf


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
function plot_ad_detections(T,out_folder,surr_time)

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
figure
set(gcf,'position',[1 1 1400 200*length(ad_rows)])
tiledlayout(length(ad_rows),1,'tilespacing','tight','padding','tight')
for ia = 1:length(ad_rows)
    nexttile
    row = ad_rows(ia);
    ad_time = currT.OnTime(row);
    ad_ch = currT.Channels{row};
    ieeg_name = currT.FileName{row};
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
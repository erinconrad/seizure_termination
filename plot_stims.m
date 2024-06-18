function plot_stims(chLabels,chs_of_interest,data_times,data_stream,...
    time_stamps,fs)


nchs = length(chs_of_interest);
figure
tiledlayout(4,1,'tilespacing','tight','padding','tight');
all_ax = nan(4,1);
for j =1:nchs
    ch = strcmp(chLabels,chs_of_interest{j});
    all_ax(j) = nexttile;
    plot(data_times,data_stream(:,ch))
    xlim([0 max(data_times)])

    if j == 4
        xlabel('time (s)')
    else
        xticklabels([])
    end
    
    ylabel('Raw signal')
    hold on
    %{
    if j == 1 || j == 2
        plot([stim_begun(j)/fs stim_begun(j)/fs],ylim,'g--')
        plot([stim_ended(j)/fs stim_ended(j)/fs],ylim,'r--')

    end

    if strcmp(ad_ch,chs_of_interest{j})
        plot([ad_begun/fs ad_begun/fs],ylim,'g--')
        plot([ad_ended/fs ad_ended/fs],ylim,'r--')
    end
    %}
end

for i = 1:size(time_stamps,1)
    ch = time_stamps{i,2};
    j = strcmp(ch,chs_of_interest);
    if sum(j) == 0, error('what'); end
    axes(all_ax(j))

    switch time_stamps{i,3}
        case 'Stim start'
            which_color = 'b';
        case 'Stim end'
            which_color = 'k';
        case 'ADs start'
            which_color = 'g';
        case 'ADs end'
            which_color = 'r';
    end
    plot([time_stamps{i,1}/fs time_stamps{i,1}/fs],ylim,'--',...
        'color',which_color,'linewidth',2)
end

end
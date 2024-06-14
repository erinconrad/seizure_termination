function plot_stims(chLabels,chs_of_interest,data_times,data_stream,...
    stim_begun,stim_ended,ad_begun,ad_ended,ad_ch,fs)
nchs = length(chs_of_interest);
figure
tiledlayout(4,1,'tilespacing','tight','padding','tight');
for j =1:nchs
    ch = strcmp(chLabels,chs_of_interest{j});
    nexttile
    plot(data_times,data_stream(:,ch))
    xlim([0 max(data_times)])

    if j == 4
        xlabel('time (s)')
    else
        xticklabels([])
    end
    
    ylabel('Raw signal')
    hold on
    if j == 1 || j == 2
        plot([stim_begun(j)/fs stim_begun(j)/fs],ylim,'g--')
        plot([stim_ended(j)/fs stim_ended(j)/fs],ylim,'r--')

    end

    if strcmp(ad_ch,chs_of_interest{j})
        plot([ad_begun/fs ad_begun/fs],ylim,'g--')
        plot([ad_ended/fs ad_ended/fs],ylim,'r--')
    end
end

end
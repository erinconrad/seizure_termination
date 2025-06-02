function aT = turn_anns_to_table(anns)

% Get annotations
n_layers = length(anns);
all_anns = {};

for l = 1:n_layers
    curr_ann_count = 0;
    while 1 % ieeg only lets you pull 250 annotations at a time
        if curr_ann_count == 0
            a=session.data.annLayer(l).getEvents(0);
        else
            a=session.data.annLayer(l).getNextEvents(a(curr_ann_count));
        end
        curr_ann_count = length(a);
        for k = 1:length(a)
            all_anns(end+1,:) = ...
                {l, a(k).start/(1e6), a(k).stop/(1e6), a(k).type, a(k).description};
        end
        if isempty(a), break; end
    end
end
aT = cell2table(all_anns,'VariableNames',{'Layer_num','Start','Stop','Type','Description'});


end
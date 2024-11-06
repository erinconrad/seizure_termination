function clean_labels = decompose_labels(chLabels)

clean_labels = cell(length(chLabels),1);

for ich = 1:length(chLabels)
    if ischar(chLabels)
        label = chLabels;
    else
        label = chLabels{ich};
    end

    %% if it's a string, convert it to a char
    if strcmp(class(label),'string')
        label = convertStringsToChars(label);
    end
    
    %% Remove leading zero
    % get the non numerical portion
    label_num_idx = regexp(label,'\d');
    if ~isempty(label_num_idx)

        label_non_num = label(1:label_num_idx-1);

        label_num = label(label_num_idx:end);

        % Remove leading zero
        if strcmp(label_num(1),'0')
            label_num(1) = [];
        end

        label = [label_non_num,label_num];
    end
    clean_labels{ich} = label;
end

end
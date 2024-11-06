function closed_ann = find_relay_anns(aT,curr_times)

% limit table to only include annotations in the specific time range
T = aT(aT.Start >= curr_times(1) & aT.Start <= curr_times(2),:);

% Look for closed relay or start stimulation
ann_text = T.Type;
closed_relay = contains(ann_text,'Closed relay') | contains(ann_text,'Start Stimulation');

if sum(closed_relay) == 0
    closed_ann = [];
    return
end

if sum(closed_relay) > 1
    error('multiple closed relays?')
end

type = ann_text{closed_relay};

% Otherwise, figure out stim channels
% See if it's a close relay
if contains(type,'Closed relay to') || contains(type, 'Start Stimulation')
    
 
    if contains(type,'Closed relay to')
    
        % find the electrodes
        C = strsplit(type);
        
        % fix for surprising text
        
        if length(C) == 6 && strcmp(C{5},'and')
            % expected order
            elec1_cell = C(end-2);
            elec2_cell = C(end);
        elseif length(C) == 8 && strcmp(C{6},'and') && ...
                (strcmp(C{4},'L') || strcmp(C{4},'R'))
            % split up L/R and rest of electrode name
            elec1_cell = {[C{4},C{5}]};
            elec2_cell = {[C{7},C{8}]};
        elseif length(C) == 8 && strcmp(C{6},'and') && ...
            (strcmp(C{4}(1),'L') || strcmp(C{4}(1),'R')) % format is "RA" "1"
            % combine RA and number
            elec1_cell = {[C{4},C{5}]};
            elec2_cell = {[C{7},C{8}]};
        elseif length(C) == 10 && strcmp(C{7},'and') && ...
                (strcmp(C{4},'L') || strcmp(C{4},'R'))
            % split up L/R and rest of electrode name
            elec1_cell = {[C{4},C{5},C{6}]};
            elec2_cell = {[C{8},C{9},C{10}]};
        else
            error('Surprising closed relay text');
            
        end
    elseif contains(type, 'Start Stimulation')

        C = strsplit(type);

        if strcmp(C{3},'from') && strcmp(C{5},'to')
            elec1_cell = C(4);
            elec2_cell = C(6);
        else

            error('surprising start stimulation text')

        end

    else
        error('what')
    end
    
    
    elec1 = elec1_cell{1};
    elec2 = elec2_cell{1};
    
   
    
end


closed_ann = sprintf('%s-%s',elec1,elec2);
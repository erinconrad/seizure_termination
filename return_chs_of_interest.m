function [chs_of_interest,stim_ch] = return_chs_of_interest(stim_pair)

C = strsplit(stim_pair,'-');
stim_ch1 = C{1};
stim_ch2 = C{2};

% Get nuerical portions
num1s = (regexp(stim_ch1, '\d+', 'match'));
num2s = (regexp(stim_ch2, '\d+', 'match'));
num1 = str2num(num1s{1});
num2 = str2num(num2s{1});
letter1 = regexp(stim_ch1,'\D+','match');
letter2 = regexp(stim_ch2,'\D+','match');
letter1 = letter1{1};letter2 = letter2{1};
assert(strcmp(letter1,letter2))


% Get lower and upper channels
num_lower = min([num1,num2])-1;
num_higher = max([num1,num2])+1;
if num_lower < 10
    str_num_lower = sprintf('0%d',num_lower);
else
    str_num_lower = sprintf('%d',num_lower);
end
if num_higher < 10
    str_num_higher = sprintf('0%d',num_higher);
else
    str_num_higher = sprintf('%d',num_higher);
end
lower_ch = [letter1,str_num_lower];
higher_ch = [letter1,str_num_higher];

chs_of_interest = {stim_ch1,stim_ch2,lower_ch,higher_ch};
stim_ch = [1 1 0 0];

end
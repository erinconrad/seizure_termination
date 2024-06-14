function is_stim = detect_stimulation(data, fs, freq)
    % Function to detect stimulation artifact (50 Hz signal)
    % Use a bandpass filter around the stimulation frequency
    [b, a] = butter(2, [freq-5 freq+5]/(fs/2), 'bandpass');
    filtered_data = filter(b, a, data);

    [c,d] = butter(2,[5 30]/(fs/2),'bandpass');
    filtered2 = filter(c,d,data);
    power = sum(filtered_data.^2);
    power2 = sum(filtered2.^2);
    rel_power = power/power2;
    threshold = 0.1; % Adjust this threshold as needed
    is_stim = rel_power > threshold;
end
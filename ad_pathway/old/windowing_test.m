fs = 1000;

windowsz = 300; % in samples
nsamples = 30000;

% sample data with 60 Hz noise added
t = (0:(nsamples-1))'/fs;
Xtrue = randn(nsamples,1);
Xnoisy = Xtrue + 10*sin(t*60*2*pi);

%plot(t,X);

%% notch filter each window separately

[b,a] = butter(4, [58 62]/(fs/2), 'stop');

Xfilter_win = zeros(nsamples,1);
for i = 1:(nsamples/windowsz)
    Xfilter_win(((i-1)*windowsz+1):(i*windowsz)) = ...
        filter(b,a,Xnoisy(((i-1)*windowsz+1):(i*windowsz)));
end

Xfilter_state = zeros(nsamples,1);
z = zeros(max(length(a),length(b))-1,1); % filter state
for i = 1:(nsamples/windowsz) 
    [Xfilter_state(((i-1)*windowsz+1):(i*windowsz)),z] = ...
        filter(b,a,Xnoisy(((i-1)*windowsz+1):(i*windowsz)),z);
end



Xfilter = filter(b,a,Xnoisy);

if 0
    plot(t,Xtrue)
    hold on
    plot(t,Xfilter_win)
    plot(t,Xfilter_state)
    legend({'true denoised','windowed notched','saved-state notched'})
else
    plot(t, (Xtrue - Xfilter_win).^2)
    hold on
    plot(t, (Xtrue - Xfilter_state).^2)
end

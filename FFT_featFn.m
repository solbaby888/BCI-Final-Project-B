function feat = FFT_featFn(signal, fs, window, overlap);

for i = 1:length(signal);
    signal_detrend = detrend(signal{i});
    
    [s w t] = spectrogram(signal_detrend, window*fs, overlap*fs, (1:60), fs, 'yaxis');
    FFT_1_60    = mean(abs(s));

    [s w t] = spectrogram(signal_detrend, window*fs, overlap*fs, (60:100), fs, 'yaxis');
    FFT_60_100   = mean(abs(s));

    [s w t] = spectrogram(signal_detrend, window*fs, overlap*fs, (100:300), fs, 'yaxis');
    FFT_100_300  = mean(abs(s));

    feat{:,i} = [FFT_1_60' FFT_60_100' FFT_100_300'];
end

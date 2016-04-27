function feat = FFT_featFn(signal, fs, window, overlap);

for i = 1:length(signal);
    signal_detrend = detrend(signal{i});
    
    [s w t] = spectrogram(signal_detrend, window*fs, overlap*fs, (5:15), fs, 'yaxis');
    FFT_5_15    = mean(abs(s));

    [s w t] = spectrogram(signal_detrend, window*fs, overlap*fs, (20:25), fs, 'yaxis');
    FFT_20_25   = mean(abs(s));

    [s w t] = spectrogram(signal_detrend, window*fs, overlap*fs, (75:115), fs, 'yaxis');
    FFT_75_115  = mean(abs(s));

    [s w t] = spectrogram(signal_detrend, window*fs, overlap*fs, (125:160), fs, 'yaxis');
    FFT_125_160 = mean(abs(s));

    [s w t] = spectrogram(signal_detrend, window*fs, overlap*fs, (160:175), fs, 'yaxis');
    FFT_160_175 = mean(abs(s));

    feat{:,i} = [FFT_5_15' FFT_20_25' FFT_75_115' FFT_125_160' FFT_160_175'];
end

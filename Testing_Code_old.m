function [pred_upsample] = Testing_Code(ECoG_Sub_Chan, fs_ECOG, no_of_channels_ECOG, weights)

%% Apply 60 Hz filter and Bandpass filter
%{
    Will eventually have to save as a function to run on all channels
%}
% 60 Hz Notch Filter - ver 2
d = designfilt('bandstopiir','FilterOrder',2, ...
           'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
           'DesignMethod','butter','SampleRate',fs_ECOG);
for i = 1:no_of_channels_ECOG;
    ECoG_Sub_Chan_filt{i} = filtfilt(d, ECoG_Sub_Chan{i});  
end;

%% Checkpt 1
%{ 
    window of 100 ms; overlap of 50 ms

    Some comments: I accidently decided to calculate some other features.
    We can remove them for now for the first checkpoint. Stupid extra work.

%}

% ====== Function handle to calculate number of windows ================= %
    NumWins = @(xLen, fs, winLen, winDisp) 1+floor((xLen-fs*winLen)/(winDisp*fs));  
    % winLen and winDisp are in secs, xLen are in samples
% ======================================================================= %

% Number of Windows for Subj
window_time = 100*10^-3;                                 % (secs)
overlap     = 50*10^-3;                                  % (secs)
windowLen   = window_time * fs_ECOG;                     % number of pts per window 
L           = length(ECoG_Sub_Chan_filt{1});
NoW         = NumWins(L, fs_ECOG, window_time, overlap); 

% Calculate frequency domain avg power features over windows
%# Calculates spectrogram and respective averages for assignment-specified frequencies
Freq_domain_Feat     = FFT_featFn(ECoG_Sub_Chan_filt, fs_ECOG, window_time, overlap); % windows X 5 Matrix

% Calculate time domain avg voltage 
kernnel = repmat(1/(window_time*fs_ECOG), 1, window_time*fs_ECOG);

for i=1:no_of_channels_ECOG;
    time_avg_volt_temp  = conv(ECoG_Sub_Chan_filt{1,i}, kernnel, 'valid');
    time_avg_volt(:, i) = time_avg_volt_temp(1:overlap*fs_ECOG:end);  
end

% Feature Matrix
for i=1:no_of_channels_ECOG
    Feat_Mat{i} = [time_avg_volt(:,i) Freq_domain_Feat{1,i}];    % time window X no. of feature
end;

% Linear Regression Prediction
numoffeat       = 6;
numofprev_win   = 3;
n_of_R          = NoW - numofprev_win;                              % number of windows for regression
p_of_R          = no_of_channels_ECOG * numoffeat * numofprev_win;            
R_mat           = zeros(n_of_R, p_of_R);                            % six features per window
curr_pt         = 3;

%# Creating R Matrix
for i = 1:n_of_R;
        curr_pt = 1 + curr_pt;
        for j = 1:no_of_channels_ECOG;
            R_idx1 = (j-1)* numoffeat * numofprev_win + 1;
            R_idx2 = R_idx1 + numoffeat * numofprev_win - 1;
            R_mat(i, R_idx1:R_idx2) = reshape(Feat_Mat{j}(curr_pt - 3:curr_pt - 1, :)', [1, 18]);
        end;
end; 

%## Adding the first columns of ones
R_ones     = ones(length(R_mat),1);
R_matrix      = [R_ones R_mat];

%# Prediction
 for i = 1:5
    pred{i}      = R_matrix*weights.f{i};
 end;

% Spline Interpolation
points_excluded = windowLen + overlap*fs_ECOG;
x  = overlap*fs_ECOG:overlap*fs_ECOG:L-points_excluded-overlap*fs_ECOG;
xx = 1:L-points_excluded;

pred_upsample = [];
for i = 1:5
    pred_upsample(:,i) = spline(x,pred{i}, xx);
end
pred_upsample = [zeros(points_excluded, 5); pred_upsample];







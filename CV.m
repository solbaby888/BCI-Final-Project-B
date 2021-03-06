    % Checkpoint 1
function [f, Correlation] = Training_Code(ECoG_Sub_Chan, fs_ECOG, no_of_channels_ECOG, Glovedata);
%{ 
    % Working only with Subject 1 but with all the channels of subject 1
    
    Table of Contents (By Section):
    % Global Vars
    % Apply 60 Hz filter
    % Checkpt 2
%}

%% Global Variables
global fs_ECOG;
global Glovedata;     
global overlap;
global windowLen;
global finger;
global L;
global subj;

%% Apply 60 Hz filter and Bandpass filter
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

% FEATURES A: Calculate frequency domain avg power features over windows
%# Calculates spectrogram and respective averages for assignment-specified frequencies
Freq_domain_Feat     = FFT_featFn_old(ECoG_Sub_Chan_filt, fs_ECOG, window_time, overlap); % windows X 5 Matrix

% FEATURES B: Calculate time domain avg voltage 
kernnel = repmat(1/(window_time*fs_ECOG), 1, window_time*fs_ECOG);
for i=1:no_of_channels_ECOG;
    time_avg_volt_temp  = conv(ECoG_Sub_Chan_filt{1,i}, kernnel, 'valid');
    time_avg_volt(:, i) = time_avg_volt_temp(1:overlap*fs_ECOG:end);  
end

% Feature Matrix
for i=1:no_of_channels_ECOG
    Feat_Mat{i} = [time_avg_volt(:,i) Freq_domain_Feat{1,i}];    % time window X no. of feature
end

% Downsampling dataglove traces                 
%# Cell array glove positions (1:5)                   
for i = 1:5;
    Glovedata_ds(:,i)  =  decimate(Glovedata(:,i), overlap*fs_ECOG); 
end;

% Linear Regression Prediction
numoffeat       = size(Feat_Mat{1},2);
numofprev_win   = 3;
n_of_R          = NoW - numofprev_win;                       % number of windows for regression
p_of_R          = no_of_channels_ECOG * numoffeat * numofprev_win;            
R_mat           = zeros(n_of_R, p_of_R);                     % six features per window
curr_pt         = 3;

for i = 1:n_of_R;
        curr_pt = 1 + curr_pt;
        for j = 1:no_of_channels_ECOG;
            R_idx1 = (j-1)* numoffeat * numofprev_win + 1;
            R_idx2 = R_idx1 + numoffeat * numofprev_win - 1;
            R_mat(i, R_idx1:R_idx2) = reshape(Feat_Mat{j}(curr_pt - 3:curr_pt - 1, :)', [1, numofprev_win*numoffeat]);
        end;
end; 

% Adding the first columns of ones


% Optimizing features

inmodel = cell(1,5);
history = cell(1,5);
f       = cell(1,5);
for finger = 1:5;
%     [inmodel{finger}, history{finger}] = sequentialfs(@optimize_function,R_mat,Glovedata_ds(5:end,finger), 'nfeatures', 50);
    [~, feature] = pca(R_mat);
    R_ones_train = ones(length(feature),1);
%     R_mat_finger = R_mat(:, inmodel{finger});
%     R_mat_finger = R_mat(:, inmodel{finger});
%     R_ones_train = ones(length(R_mat_finger),1);
%     TRAIN        = [R_ones_train R_mat_finger];
    TRAIN        = [R_ones_train feature];
    f{finger}    = mldivide(TRAIN' * TRAIN, TRAIN' * Glovedata_ds(5:end,finger));
    pred           = TRAIN*f{finger}; 
    
%     spline
    points_excluded = windowLen + overlap*fs_ECOG;
    x               = overlap*fs_ECOG:overlap*fs_ECOG:L-points_excluded-overlap*fs_ECOG;
    xx              = 1:L-points_excluded;
    pred_upsample   = spline(x, pred, xx);
    pred_upsample   = [zeros(points_excluded, 1); pred_upsample'];
   
    Correlation(finger)     = corr(Glovedata(:,finger), pred_upsample);
end;
% Save
save(strcat('inmodel', '_sub', num2str(subj)), 'inmodel');
save(strcat('history', '_sub', num2str(subj)), 'history');
save(strcat('weights', '_sub', num2str(subj)), 'f');



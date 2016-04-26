function criterion = my_function(XTRAIN, YTRAIN, XTEST, YTEST) %(XTRAIN,y_train,X_test,y_test)

% XTRAIN
% Linear Regression Prediction

% split up xtrain into component ecog feature matrices. Use individual
% matrices to create R matrix? 
% 

numoffeat       = size(XTRAIN{1},2);
numofprev_win   = 3;
n_of_R          = size(XTRAIN{1},1) - numofprev_win;                       % number of windows for regression
p_of_R          = size(XTRAIN) * numoffeat * numofprev_win;                % SO CONFUSED NEED TO FIX         
R_mat           = zeros(n_of_R, p_of_R);                                   % six features per window
curr_pt         = 3;

for i = 1:n_of_R;
        curr_pt = 1 + curr_pt;
        for j = 1:no_of_channels_ECOG;
            R_idx1 = (j-1)* numoffeat * numofprev_win + 1;
            R_idx2 = R_idx1 + numoffeat * numofprev_win - 1;
            R_mat(i, R_idx1:R_idx2) = reshape(XTRAIN{j}(curr_pt - 3:curr_pt - 1, :)', [1, numofprev_win*numoffeat]);
        end;
end; 

% Adding the first columns of ones
R_ones     = ones(length(R_mat),1);
R_matrix   = [R_ones R_mat];

% Weights and prediction for each finger

f{i}         = mldivide(R_matrix'*R_matrix, R_matrix' * Glovedata_ds{i}(5:end));
% pred{i}      = R_matrix*f{i};   end;

% Spline interpolation
%# the interpolation isnt working yet..needs work
% points_excluded = windowLen + overlap*fs_ECOG;
% x  = overlap*fs_ECOG:overlap*fs_ECOG:L-points_excluded-overlap*fs_ECOG;
% xx = 1:L-points_excluded;

% pred_upsample = [];
% for i = 1:5
%     pred_upsample(:,i) = spline(x, pred{i}, xx);
% end
% pred_upsample = [zeros(points_excluded, 5); pred_upsample];
% for i = 1:5;
%    Glovedata_mat(:,i) = Glovedata{i};
% end;

% XTEST CALCULATION

save(strcat('weight', '_CV'), 'f');

Training_Correlation = corr(Glovedata_mat, round(pred_upsample));

criterion = Training_Correlation;


end
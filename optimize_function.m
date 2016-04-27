function criterion = optimize_function(XTRAIN, YTRAIN, XTEST, YTEST) %(XTRAIN,y_train,X_test,y_test)
 global overlap;
 global windowLen;
 global finger;
 global L;
 global fs_ECOG;
 global Glovedata;
 global subj; 
 
% XTRAIN
% Linear Regression Prediction
% fun takes different features of R matrix. No of rows remains consistent.
% Weights and prediction for each finger

f        = mldivide(XTRAIN' * XTRAIN, XTRAIN' * YTRAIN);
pred     = XTEST*f; 

%Spline interpolation
%# the interpolation isnt working yet..needs work
% points_excluded = windowLen + overlap*fs_ECOG;
% x               = overlap*fs_ECOG:overlap*fs_ECOG:L-points_excluded-overlap*fs_ECOG;
% xx              = 1:L-points_excluded;
% 
% 
% pred_upsample = spline(x, pred, xx);
% pred_upsample = [zeros(points_excluded, 1); pred_upsample];

% XTEST CALCULATION
Training_Err = -corr(YTEST, pred);
criterion    = Training_Err;

end
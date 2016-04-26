function criterion = my_function(XTRAIN, YTRAIN, XTEST, YTEST) %(XTRAIN,y_train,X_test,y_test)

% XTRAIN
% Linear Regression Prediction

% fun takes different features of R matrix. No of rows resmain consistent.


% Weights and prediction for each finger

f        = mldivide(XTRAIN' * XTRAIN, XTRAIN' * YTRAIN);
pred      = XTEST*f; 

%Spline interpolation
%# the interpolation isnt working yet..needs work
points_excluded = windowLen + overlap*fs_ECOG;
x  = overlap*fs_ECOG:overlap*fs_ECOG:L-points_excluded-overlap*fs_ECOG;
xx = 1:L-points_excluded;

pred_upsample = [];
for i = 1:5
    pred_upsample(:,i) = spline(x, pred{i}, xx);
end
pred_upsample = [zeros(points_excluded, 5); pred_upsample];

% XTEST CALCULATION
% 
% save(strcat('weight', '_CV'), 'f');

% need to figure out how to pass specific glovedata for each finger
% iteration & subject.
Training_Correlation = corr(Glovedata{i},pred_upsample);

criterion = Training_Correlation;


end
%% Training
%{ 
    % Working only with all Subjects
    
    Table of Contents (By Section):
    % Import Data
    % Training on subj 1 2 3
%}
close all;
clear all;
clc;

% Import Data
% ==Login=================================================================%
% Karan's login
% ieeglogin = 'kar_ieeglogin.bin';
% username = 'karanm';
% Thuy's login
username = 'solbaby888';
ieeglogin = 'sol_ieeglogin.bin';
%=========================================================================%

% Subject 1 2 3
for subj = 1:3;
    
    if subj == 1;
        filename = 'I521_A0012';
    elseif subj == 2;
        filename = 'I521_A0013';
    elseif subj == 3;
        filename = 'I521_A0014';
    end;
    
    sesh_sub_1 = IEEGSession(strcat(filename, '_D001'), username, ieeglogin);
    sesh_sub_2 = IEEGSession(strcat(filename, '_D002'), username, ieeglogin);
    sesh_sub_3 = IEEGSession(strcat(filename, '_D003'), username, ieeglogin);


    % Parameters for EcoG data
    fs_ECOG             = sesh_sub_1.data(1).sampleRate;            % Sampling Rate
    no_of_channels_ECOG = size(sesh_sub_1.data(1).rawChannels, 2);  % Number of Channels 
  
    % ECoG data
    nr          = ceil((sesh_sub_1.data(1).rawChannels(1).get_tsdetails.getEndTime)/...
                1e6*sesh_sub_1.data(1).sampleRate);
    ECoG_Sub_Chan       = cell(1, no_of_channels_ECOG);
    for i = 1:no_of_channels_ECOG;
        ECoG_Sub_Chan{:, i} = smooth(sesh_sub_1.data(1).getvalues(1:nr, i)); 
    end;
    
    % Glove data       
    nr        = ceil((sesh_sub_2.data(1).rawChannels(1).get_tsdetails.getEndTime)/...
                            1e6*sesh_sub_2.data(1).sampleRate);
    Glovedata = {};
    for i = 1:5;
        Glovedata{i} = sesh_sub_2.data(1).getvalues(1:nr, i);
    end;

    % Parameters of Glove data
    fs_dataglove        = sesh_sub_2.data(1).sampleRate;
    
    % TRAIN
    [f, ~, Training_Correlation] = Training_Code(ECoG_Sub_Chan, fs_ECOG, no_of_channels_ECOG, Glovedata);
    
    % SAVE
    % Save weights;
    save(strcat('weights', '_sub', num2str(subj)), 'f');
    % Save training correlation
    save(strcat('training', '_sub', num2str(subj)), 'Training_Correlation');
    % Save predictions
    % save(strcat('pred', '_sub', num2str(subj)), 'pred_upsample');

   
end;  


%% Cross Validation of Subj 1
% clear;clc;
% %% Import Data
% % ==Login=================================================================%
% % Karan's login
% % ieeglogin = 'kar_ieeglogin.bin';
% % username = 'karanm';
% % Thuy's login
% username = 'solbaby888';
% ieeglogin = 'sol_ieeglogin.bin';
% %=========================================================================%
% 
% % Subject 1
% filename = 'I521_A0012';
% sesh_sub_1 = IEEGSession(strcat(filename, '_D001'), username, ieeglogin);
% sesh_sub_2 = IEEGSession(strcat(filename, '_D002'), username, ieeglogin);
% 
%  
% % Parameters for EcoG data
% fs_ECOG             = sesh_sub_1.data(1).sampleRate;            % Sampling Rate
% no_of_channels_ECOG = size(sesh_sub_1.data(1).rawChannels, 2);  % Number of Channels 
% 
% 
% % ECoG data
% nr          = ceil((sesh_sub_1.data(1).rawChannels(1).get_tsdetails.getEndTime)/...
%             1e6*sesh_sub_1.data(1).sampleRate);
% ECoG_Sub_Chan       = cell(1, no_of_channels_ECOG);
% for i = 1:no_of_channels_ECOG;
%     ECoG_Sub_Chan{:, i} = smooth(sesh_sub_1.data(1).getvalues(1:nr, i)); 
% end;
% L                   = length(ECoG_Sub_Chan{1});
% 
% % Glove data       
% nr        = ceil((sesh_sub_2.data(1).rawChannels(1).get_tsdetails.getEndTime)/...
%                         1e6*sesh_sub_2.data(1).sampleRate);
% Glovedata = {};
% for i = 1:5;
%     Glovedata{i} = sesh_sub_2.data(1).getvalues(1:nr, i);
% end;
% overlap     = 50*10^-3;
% for i = 1:5;
%     Glovedata_ds{i}  =  decimate(Glovedata{i}, overlap*fs_ECOG); 
% end;
% 
% [f, ~, Training_Correlation, R_matrix] = Training_Code(ECoG_Sub_Chan, fs_ECOG, no_of_channels_ECOG, Glovedata);
% 
% % CV
% indices = crossvalind('Kfold',length(R_matrix),10);
% Testing_Correlation = cell(10,1);
% for i = 1:10
%     test  = (indices == i); 
%     train = ~test;
%     
%     train_R = R_matrix(train==1, :);
%     glove_train = cell(1, 5);
%     idx       = 0;
%     % Train
%     for k = 1:length(train)
%         if train(k) == 1;
%             idx             = idx + 1;
%             for j = 1:5
%                 glove_train{j}(idx, 1)  = Glovedata_ds{j}(k);
%             end;
%         end;
%     end;
%     
%     for k = 1:5
%         f{k}         = mldivide(train_R'*train_R, train_R' * glove_train{i}(5:end));
%     end;
%     
%     points_excluded = windowLen + overlap*fs_ECOG;
%     x  = overlap*fs_ECOG:overlap*fs_ECOG:L-points_excluded-overlap*fs_ECOG;
%     xx = 1:L-points_excluded;
%     
%     pred_upsample = Testing_Code(ECoG_Sub_Chan, fs_ECOG, sum(test), f);
%     Glovedata_mat = zeros(length(Glovedata{1}), 5);
%     for j = 1:5;
%         Glovedata_mat(:,j) = glove_test{j};
%     end;
%     Testing_Correlation{j} = corr(Glovedata_mat(test), round(pred_upsample));
% end

%% Testing on Testing Set
% ==Login=================================================================%
% Karan's login
% ieeglogin = 'kar_ieeglogin.bin';
% username = 'karanm';
% Thuy's login
username = 'solbaby888';
ieeglogin = 'sol_ieeglogin.bin';
%=========================================================================%
%Import Weights
% weights_sub1 = load('weights_sub3');
for subj = 1:3;
    if subj == 1;
        filename = 'I521_A0012';
        weights = load('weights_sub1');
    elseif subj == 2;
        filename = 'I521_A0013';
        weights = load('weights_sub2');
    elseif subj == 3;
        filename = 'I521_A0014';
        weights = load('weights_sub3');
    end;

    sesh_sub_3 = IEEGSession(strcat(filename, '_D003'), username, ieeglogin);
    nr          = ceil((sesh_sub_3.data(1).rawChannels(1).get_tsdetails.getEndTime)/...
                1e6*sesh_sub_3.data(1).sampleRate);
            

    % Parameters for EcoG data
    fs_ECOG             = sesh_sub_3.data(1).sampleRate;            % Sampling Rate
    no_of_channels_ECOG = size(sesh_sub_3.data(1).rawChannels, 2);  % Number of Channels 

    % Testing with all channels of each subject
    no_of_channels_ECOG = size(sesh_sub_3.data(1).rawChannels, 2); 
    ECoG_Sub_Chan       = cell(1, no_of_channels_ECOG);

    for i = 1:no_of_channels_ECOG;
        ECoG_Sub_Chan{:, i} = smooth(sesh_sub_3.data(1).getvalues(1:nr, i)); 
    end;
    
    pred_upsample = Testing_Code(ECoG_Sub_Chan, fs_ECOG, no_of_channels_ECOG, weights);
	save(strcat('pred', '_sub', num2str(subj)), 'pred_upsample');
    
    % Submission Format
    predicted_dg{subj,1}=(pred_upsample);
    save('predicted_dg', 'predicted_dg');
end;

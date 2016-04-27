%% Training
%{ 
    % Working only with all Subjects
    
    Table of Contents (By Section):
    % Import Data
    % Training on subj 1 2 3
%}
% close all;
% clear all;
% clc;

% Import Data
% ==Login=================================================================%
% Karan's login
% ieeglogin = 'kar_ieeglogin.bin';
% username = 'karanm';
% Thuy's login
username = 'solbaby888';
ieeglogin = 'sol_ieeglogin.bin';
%=========================================================================%
% Save Datasets
% Subject 1 2 3
% for subj = 1:3;
%     global Glovedata;
%     if subj == 1;
%         filename = 'I521_A0012';
%     elseif subj == 2;
%         filename = 'I521_A0013';
%     elseif subj == 3;
%         filename = 'I521_A0014';
%     end;
%     
%      sesh_sub_1 = IEEGSession(strcat(filename, '_D001'), username, ieeglogin);
%      sesh_sub_2 = IEEGSession(strcat(filename, '_D002'), username, ieeglogin);
%      sesh_sub_3 = IEEGSession(strcat(filename, '_D003'), username, ieeglogin);
% 
%    no_of_channels_ECOG = size(sesh_sub_1.data(1).rawChannels, 2);  % Number of Channels 
% 
%     % ECoG data
%     nr                  = ceil((sesh_sub_1.data(1).rawChannels(1).get_tsdetails.getEndTime)/...
%                             1e6*sesh_sub_1.data(1).sampleRate);
%     ECoG_Sub_Chan       = cell(1, no_of_channels_ECOG);
%     for i = 1:no_of_channels_ECOG;
%         ECoG_Sub_Chan{:, i} = smooth(sesh_sub_1.data(1).getvalues(1:nr, i)); % SMOOTHED DATA: May want to optimize
%     end;
%     
%     % Glove data       
%     nr        = ceil((sesh_sub_2.data(1).rawChannels(1).get_tsdetails.getEndTime)/...
%                             1e6*sesh_sub_2.data(1).sampleRate);
%     Glovedata = [];
%     for i = 1:5;
%         Glovedata(:, i) = sesh_sub_2.data(1).getvalues(1:nr, i);
%     end;
%     
%     save(strcat('Sub', num2str(subj), '_Channels'), 'ECoG_Sub_Chan')
%     save(strcat('Sub', num2str(subj), '_Glovedata'), 'Glovedata')
%       
% end;  
%=========================================================================%
% Subject 1 2 3
global subj;
for subj = 1:3;
    % Global Variables
    global fs_ECOG;
    global Glovedata;     
    global overlap;
    global windowLen;
    global finger;
    global L;
    global subj;
    
    if subj == 1;
        filename = 'I521_A0012';
        load('Sub1_Channels')
        load('Sub1_Glovedata')
    elseif subj == 2;
        filename = 'I521_A0013';
        load('Sub2_Channels')
        load('Sub2_Glovedata')
    elseif subj == 3;
        filename = 'I521_A0014';
        load('Sub3_Channels')
        load('Sub3_Glovedata')
    end;
    
     sesh_sub_1 = IEEGSession(strcat(filename, '_D001'), username, ieeglogin);
     sesh_sub_2 = IEEGSession(strcat(filename, '_D002'), username, ieeglogin);
     sesh_sub_3 = IEEGSession(strcat(filename, '_D003'), username, ieeglogin);

    % Parameters of Glove data
     fs_dataglove        = sesh_sub_2.data(1).sampleRate;
    % Parameters for EcoG data
     fs_ECOG             = sesh_sub_1.data(1).sampleRate;            % Sampling Rate
     no_of_channels_ECOG = size(sesh_sub_1.data(1).rawChannels, 2);  % Number of Channels 
 
    % TRAIN
    [f, Correlation] = Training_Code(ECoG_Sub_Chan, fs_ECOG, no_of_channels_ECOG, Glovedata);     
    save(strcat('Training_correlation', '_sub', num2str(subj)), 'Correlation'); 
end;  

%%      


%% Cross Validation of Subj 1
% clear;clc;
% % Import Data
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
% % Parameters for EcoG data
% fs_ECOG             = sesh_sub_1.data(1).sampleRate;            % Sampling Rate
% no_of_channels_ECOG = size(sesh_sub_1.data(1).rawChannels, 2);  % Number of Channels 
% 
% % ECoG data
% nr          = ceil((sesh_sub_1.data(1).rawChannels(1).get_tsdetails.getEndTime)/...
%             1e6*sesh_sub_1.data(1).sampleRate);
% ECoG_Sub_Chan_train       = cell(1, no_of_channels_ECOG);
% ECoG_Sub_Chan_test        = cell(1, no_of_channels_ECOG);
% for i = 1:no_of_channels_ECOG;
%     ECoG_Sub_Chan_train{:, i} = sesh_sub_1.data(1).getvalues(1:nr*3/5, i); 
%     ECoG_Sub_Chan_test{:, i}  = sesh_sub_1.data(1).getvalues((nr*2/5 + 1):nr, i);
% end;
% L                   = length(ECoG_Sub_Chan_train{1});
% 
% % Glove data       
% nr        = ceil((sesh_sub_2.data(1).rawChannels(1).get_tsdetails.getEndTime)/...
%                         1e6*sesh_sub_2.data(1).sampleRate);
% Glovedata_train = {};
% Glovedata_test  = {};
% for i = 1:5;
%     Glovedata_train{i} = sesh_sub_2.data(1).getvalues(1:nr*3/5, i);
%     Glovedata_test{i}  = sesh_sub_2.data(1).getvalues((nr*2/5 + 1):nr, i);
% end;
% 
% [f, ~, Training_Correlation, R_matrix] = Training_Code(ECoG_Sub_Chan_train, fs_ECOG, no_of_channels_ECOG, Glovedata_train);
% save(strcat('weight', '_CV'), 'f');
% weight = load('weight_CV');
% pred_upsample                          = Testing_Code(ECoG_Sub_Chan_test, fs_ECOG, no_of_channels_ECOG, weight);
% 
% for j = 1:5;
%         Glovedata_mat(:,j) = Glovedata_test{j};
% end;
% Testing_Correlation = corr(Glovedata_mat, round(pred_upsample));
% 
% % CV
% % indices = crossvalind('Kfold',length(R_matrix),10);
% % Testing_Correlation = cell(10,1);
% % for i = 1:10
% %     test  = (indices == i); 
% %     train = ~test;
% %     
% %     train_R = R_matrix(train==1, :);
% %     glove_train = cell(1, 5);
% %     idx       = 0;
% %     % Train
% %     for k = 1:length(train)
% %         if train(k) == 1;
% %             idx             = idx + 1;
% %             for j = 1:5
% %                 glove_train{j}(idx, 1)  = Glovedata_ds{j}(k);
% %             end;
% %         end;
% %     end;
% %     
% %     for k = 1:5
% %         f{k}         = mldivide(train_R'*train_R, train_R' * glove_train{i}(5:end));
% %     end;
% %     
% %     points_excluded = windowLen + overlap*fs_ECOG;
% %     x  = overlap*fs_ECOG:overlap*fs_ECOG:L-points_excluded-overlap*fs_ECOG;
% %     xx = 1:L-points_excluded;
% %     
% %     pred_upsample = Testing_Code(ECoG_Sub_Chan, fs_ECOG, sum(test), f);
% %     Glovedata_mat = zeros(length(Glovedata{1}), 5);
% %     for j = 1:5;
% %         Glovedata_mat(:,j) = glove_test{j};
% %     end;
% %     Testing_Correlation{j} = corr(Glovedata_mat(test), round(pred_upsample));
% % end

%% Testing on Testing Set
% ==Login=================================================================%
% Karan's login
% ieeglogin = 'kar_ieeglogin.bin';
% username = 'karanm';
% Thuy's login
username = 'solbaby888';
ieeglogin = 'sol_ieeglogin.bin';
%=========================================================================%
% Save datasets
for subj = 1:3;
    if subj == 1;
        filename = 'I521_A0012';
    elseif subj == 2;
        filename = 'I521_A0013';
    elseif subj == 3;
        filename = 'I521_A0014';
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
    
    save(strcat('Sub', num2str(subj), '_Channels_Test', 'ECoG_Sub_Chan'));
end;
%=========================================================================%
% Gets Predictions
for subj = 1:3;
    %Import Weights
    if subj == 1;
        load('Sub1_Channels_Test')
        weights = load('weights_sub1');
    elseif subj == 2;
        load('Sub1_Channels_Test')
        weights = load('weights_sub2');
    elseif subj == 3;
        load('Sub1_Channels_Test')
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
    
    save(strcat('Sub', num2str(subj), '_Channels_Test'), 'ECoG_Sub_Chan');    
    pred_upsample = Testing_Code(ECoG_Sub_Chan, fs_ECOG, no_of_channels_ECOG, weights);
  
    
    save(strcat('pred', '_sub', num2str(subj)), 'pred_upsample');
   
    % Submission Format
    predicted_dg{subj,1} = (pred_upsample);
    save('predicted_dg', 'predicted_dg');
end;

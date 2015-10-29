clear all
close all;
addpath(genpath('C:\Users\havm2701\Dropbox\PhD\Matlab\MATLAB\MHA'))
global T1C T2 FLAIR 
brain_dir ='N:\User\Research\Brain_data\BRATS-2\Image_Data\';
result_path = 'N:\User\Research\Brain_data\BRATS-2_nii\Image_Data\';
if exist(result_path,'dir')==0
    mkdir(result_path)
end


%% define paths and parameters

   
Nb_Brains = 30;
brain_type_list = dir(brain_dir);
brain_type_list = brain_type_list(3:end);
 for subfoldercounter=1:length(brain_type_list)
     type = brain_type_list(subfoldercounter).name;
     brain_type_path = [brain_dir,type];
     brain_list = dir(brain_type_path);
     brain_list =brain_list(3:end);

     for brain_counter =1:length(brain_list)
          
         name = brain_list(brain_counter).name;
         brain_name_path = [brain_type_path,'\',name];
 result_path_brain = [result_path,type,'_',name];        
if exist(result_path_brain,'dir')==0
    mkdir(result_path_brain)
end
 mask_path = ['N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\',type,'\',name];
 MASK = load_mask(mask_path);
 MSDK = double(MASK);
t1 = tic();
%% load modalities
        [T1C, T2, FLAIR, truth,T1, info,flair_name] = load_modalities(brain_name_path,1);
        [height width depth] = size(T1C);
        T1C = double(T1C);
        T2 = double(T2);
        FLAIR = double(FLAIR);
        % take backup of the modalities
        truth_backup = truth;





% disp('***********************************')
% disp(['Processing brain ',type,'_',name, ' with k = ', num2str(2)])
% disp('***********************************')
% T1C_nii = make_nii(T1C);
% save_nii(T1C_nii, [result_path_brain,'\',type,'_',name,'_T1C.nii.gz']);
% T2_nii = make_nii(T2);
% save_nii(T2_nii, [result_path_brain,'\',type,'_',name,'_T2.nii.gz']);
% FLAIR_nii = make_nii(FLAIR);
% save_nii(FLAIR_nii, [result_path_brain,'\',type,'_',name,'_FLAIR.nii.gz']);
% T1_nii = make_nii(T1);
% save_nii(T1_nii, [result_path_brain,'\',type,'_',name,'_T1.nii.gz']);
% truth_nii = make_nii(truth);
% save_nii(truth_nii, [result_path_brain,'\',type,'_',name,'_truth.nii.gz']);

mask_nii = make_nii(MASK);
save_nii(mask_nii, [result_path_brain,'\',type,'_',name,'_mask.nii.gz']);




     end
 end

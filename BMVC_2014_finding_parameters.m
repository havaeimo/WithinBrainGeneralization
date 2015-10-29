clear all
close all;
addpath(genpath('C:\Users\havm2701\Dropbox\PhD\Matlab\MATLAB\MHA'))
%run('C:\Users\havm2701\Dropbox\PhD\Matlab\MATLAB\vlfeat-0.9.17\toolbox\vl_setup');
%mex -lut C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\markov_network_multilable_general3Dgraph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\graph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\LinkedBlockList.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\maxflow.cpp
%mex -lut C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\interactiveGraphcut_withdataterms.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\GCoptimization.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\graph.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\LinkedBlockList.cpp C:\Users\havm2701\Dropbox\PhD\Matlab\mex_test\gcmex-2.3.0\GCMex\maxflow.cpp
global T1C T2 FLAIR 
%classifier_name = 'Tree';
classifier_name = 'kNN';
% classifier_name = 'SVM';
spatial_feature = 1; % set to one if you want  the spatial features
k=3;
parameters =[];
if ~strcmp(classifier_name,'kNN')
    k = NaN;
end

% output directories
if spatial_feature==0
    result_sub_directory = 'Newmasks_T1c_T2_flair';
else
    result_sub_directory = 'Newmasks_6Dfeature';
end


%brain_dir = '/Volumes/Expansion Drive/User/Research/Brain_data/BRATS-2/Image_Data/';
%result_path_brain = 'N:\User\Research\Results\textfiles_oct\';
%mask_dir = '/Users/uoft/Dropbox/PhD/BratsAnalysis/CHALLENGE2013_KNN_mha/';


%result_path = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\KNN_stat_result_for_k\';
% result_path = 'J:\User\Research\Results\BRATS-2_KNN_ENDAUGUST_lowresT2\';

% 
% if exist(textfiles_save_dir,'dir')==0
%     mkdir(textfiles_save_dir)
% end
% 
% save_path1 = [textfiles_save_dir,'summery_evaluation_methods.txt'];
% % if exist(result_path_brain,'dir')==0
% %     mkdir(result_path_brain)
% % end

%% define paths and parameters
Nb_Selected = 1000;

Nb_healthy = 1000;
downsamplerate = 20;




%  brain_name_path = ['N:\User\Research\Brain_data\BRATS-2\Image_Data\',type,'\',name];
%  mask_name_path = ['N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\',type,'\',name];

 mask_dir = ['N:\User\Research\Brain_data\BRATS-2_MASKS\Image_Data\'];
 brain_dir ='N:\User\Research\Brain_data\BRATS-2\Image_Data\';
%mask_dir = 'N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\';
%alpha_vector = 0.03:0.01:0.1;
%alpha_vector =  0.003:0.005:.05;

%alpha_vector =  0.0005:0.001:.008;
alpha_vector = 0.02:0.005:0.1;
beta_vector =  0.0045 : 0.0003 : 0.006+0.0005;
% for a =1:10; 
%     for b =1:20


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
          if type == 'LG'
%              error('finished')
               continue;
          end
t1 = tic();
%% load modalities
        [T1C, T2, FLAIR, truth,T1,info,flair_name]  = load_modalities(brain_name_path,1);
        [height width depth] = size(T1C);
        T1C = double(T1C);
        T2 = double(T2);
        FLAIR = double(FLAIR);
        % take backup of the modalities
        truth_backup = truth;

%% reduce the resolution of T2 (check for Maxime)
%        T2 = rescalingT2(T2);
        
%% creating the mask and ROI
   % MASK = make_mask(T1C, T2, FLAIR);
% or load the mask from a path
%mask_name_path = 'N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\HG\0026';
 %mask_name_path = ['N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\',type,'\',name];
 mask_name_path = [mask_dir,type,'\',name];
 MASK = load_mask(mask_name_path);
 backup_MASK = MASK;
% you can Edit the mask
     %MASK = edit_mask(T1C, T2, FLAIR, MASK);
 
%[ymin,ymax,xmin,xmax,zmin,zmax] = make_regtangle(T1C);  
% save ROI (regtangle) coordinate points
 


%% Detect the borders of the brain

[ymin,ymax,xmin,xmax,zmin,zmax] = find_boundingbox_borders(T1C);


%%% load reg coordinates from the saved location  => this method is
%%% obseleate in the new version. In the new version we detect the borders
%%% of the brain by find_boundingbox_borders function 


% roi_save_dir = ['N:\User\Research\Brain_data\BRATS-2_roi_coordinates\Image_Data\',type,'\',name];
%  roi_save_path = [roi_save_dir,'\ROI_',type,'_',name,'.mat'] ; 
% load(roi_save_path)
%% 
% %  if exist(roi_save_dir,'dir')==0
%     mkdir(roi_save_dir)
%  end
%    save(roi_save_path,'ymin','ymax','xmin','xmax','zmin','zmax');

 
% mask_save_path = [mask_save_dir,'\MASK_',type,'_',name,'.mha'] ;     
% MASK = uint16(MASK);
% writemetaimagefile(mask_save_path, MASK, info.PixelDimensions,info.Offset);     
% MASK = double(MASK);  
disp('***********************************')
disp(['Processing brain ',type,'_',name, ' with k = ', num2str(k)])
disp('***********************************')






%%
T1C = T1C(xmin:xmax , ymin:ymax , zmin:zmax);
T2 = T2(xmin:xmax , ymin:ymax , zmin:zmax);
FLAIR = FLAIR(xmin:xmax , ymin:ymax , zmin:zmax);

MASK = MASK(xmin:xmax , ymin:ymax , zmin:zmax);
truth = truth(xmin:xmax , ymin:ymax , zmin:zmax);
truth = double(truth); 
%%  make space and selected points space( MADE CHANGES FOR MAKING TEXT FILES >> FIX LATER)
  
[space,mask_idx,truth_idx] = make_space(T1C, FLAIR, T2,truth, MASK);
 space2 = space;
         background = find(sum(space(1:3,:))< 0.001); 
         space(:,background) = [];
         mask_idx(background) = [];
         truth_idx(background) = [];


% Remove the spatial features
if spatial_feature ==0
    cordinates = space(4:6,:);
    space = space(1:3,:);
end


[selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );


%  make_txt_file(space,selected_space,truth_idx,mask_idx,type, name);

  
 [h,w,d] = size(T1C);
% 
%% KNN search classifier
if  strcmp(classifier_name,'kNN')
t2 = tic();
kkk=k;
%[IDX,D] = knnsearch(selected_space',space','K', kkk ,'NSMethod','kdtree' ,'Distance','euclidean');
[IDX,D] = KNN_search_rbf(selected_space', space', 3,0.5);
knn_time = toc(t2);
%[ind,knndist] = kdtreeKNN(selected_space,mask_idx,space,k);
%IDX = ind';
%D = knndist';

% assign labels
if spatial_feature ==0
    space = [space;cordinates];
end
% segmented_volume = assign_labels(IDX,D,mask_idx, space, h, w, d);
% 
% segmented_volume = medfilt3(segmented_volume, [3,3,3]);
% % compute statistics 
% segmented_brain = zeros(height,width,depth);
% segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
% 
% %  segmented_brain(segmented_brain == 3)=4;
% %  segmented_brain(segmented_brain == 1)=3;
%  
%  mha_result_name = strrep(flair_name , 'KNN',['Seg_',type,'_',name]);
% mha_result_path = [mha_result_name];
% segmented_brain = uint16(segmented_brain);
% writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset);

% make data-term_KNN

matrix = make_dataterm_matrix(IDX,mask_idx,space,T1C);
matrix = (matrix/k);
%energy = -log(matrix);
energy_base = 1 - matrix;
%% 
elseif  strcmp(classifier_name,'Tree')

    % Discriminative classifier
 t1 = tic();
    mask_idx(mask_idx==10)=0;
    ens = fitensemble(selected_space',mask_idx','LPBoost',100,'tree');
    [output,score]=predict(ens,space');
tree_time =toc(t1);
display(['tree time = ',num2str(tree_time)])
    if spatial_feature ==0
        space = [space;cordinates];
    end

    segmented_volume = assign_labels_classifier(output,mask_idx, space, h, w, d);
    matrix = make_dataterm_matrix_classifier(score,mask_idx,space,T1C);
    energy_base = 1 - matrix;
%%
elseif  strcmp(classifier_name,'SVM')
% %%%%%%%%KERNAL SVM%%%%%%%%%%%
    %inputdir = 'N:\User\Research\Results\MICCAI_2014_Feb\FEB_results\SVM\Probabilities\';
   %inputdir = 'N:\User\Research\Results\SVM_OUTPUTS_fromPhilippe\libsvm_probabilities\training_test\'
    inputdir = 'N:\User\Research\Results\SVM_OUTPUTS_fromPhilippe\libsvm_probabilities\BMVC\libsvm_training_6dim\libsvm_results\';
   %inputdir = 'N:\User\Research\Results\MICCAI_2014_Feb\FEB_results\SVM\libsvm_training_3dim\libsvm_results\';
    %inputdir = 'N:\User\Research\Results\SVM_OUTPUTS_fromPhilippe\libsvm_probabilities\challenge\libsmv_prob_6dim\'
    inputname = [type,'_',name,'_libsvm_output.txt'] ;
    inputfile = [inputdir,inputname];
    [segmented_volume posterior_matrix] = visualize_svm_outputs_2(inputfile,info,1,backup_MASK); % 1 is to get the probabilities



    posterior_matrix = posterior_matrix(: , xmin:xmax , ymin:ymax , zmin:zmax);
    segmented_volume = segmented_volume(xmin:xmax , ymin:ymax , zmin:zmax);
    energy_base = 1 - posterior_matrix;
    %save the energy matrix in .mat file format to be used later
    energy_mat_dir = [inputdir, 'matfiles\'];
    energy_mat_path = [energy_mat_dir,'energy_term_',type,'_',name,'.mat'];

    if exist(energy_mat_dir,'dir')==0
        mkdir(energy_mat_dir)
    end
   
    %save energy term to mat file
   save(energy_mat_path,'energy_base');
   % load energy term from mat file
%    var = load(energy_mat_path,'energy_base');
%    energy_base = var.energy_base;
end

%% APPLY graphcut BOYKOV JOLLY with dataterms
for a =1:length(alpha_vector); 
   %for b =1:length(beta_vector) 
       alpha = alpha_vector(a);
    %   beta = beta_vector(b);
     beta = 0.0053 ;
%        if beta == 0.004
%             continue
%         end
%         alpha = 0.01
%         beta = 0.005
disp('***********************************')
disp(['Processing brain ',type,'_',name, ' with alpha: ', num2str(alpha),'beta: ',num2str(beta)])
disp('***********************************')       
        
 result_dir =['N:\User\Research\Results\BMVC_2014_find_parameters\',classifier_name,'_rbf\',result_sub_directory,'\'];
%main_result_dir = [result_dir,'a',num2str(alpha),'b',num2str(beta),'\'];
% save_path3 = [result_dir,'parameter_stats.txt'];                  
% textfiles_save_dir =[main_result_dir,'textresult\'];

%alpha = 0.001;
 energy = alpha * energy_base;
 %beta = .0051;
 t3 = tic();
lf = interactiveGraphcut_withdataterms(space2,double(MASK),energy,beta);
graphcut_time = toc(t3)
[h,w,d] = size(T1C);
  lf = reshape(lf,h, w, d);
   segmented_volume = lf ;
   % apply median filter
segmented_volume = medfilt3(segmented_volume, [3,3,3]);
% compute statistics 
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;

 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;
 
%  mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = [mha_result_name];
% segmented_brain = uint16(segmented_brain);
% writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset);
 
[dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(segmented_brain,truth_backup);
 stats = [precision;recall;dice;jaccard];
% avg = sum(stats(:))/16;
% sum = 0;
% sub = 0;
% for l=1:length(stats(:))
%     st = stats(:);
%     val = stats(a);
%     if ~isnan(val)
%         sum = sum+val;
%         sub = sub+1;
%     end
% end
% avg = sum/(length(st)-sub);


%%%%%%%%%%%% save results
% mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = [boykov_result_save_dir,mha_result_name];
% segmented_brain = uint16(segmented_brain);
% writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset); 
%%%% save the matfile for statistics

boykov_denoised_result_save_dir =  [result_dir,'boykov\Denoised\'];
if exist(boykov_denoised_result_save_dir,'dir')==0
    mkdir(boykov_denoised_result_save_dir)
end
evaluation_matrix_savedir = [boykov_denoised_result_save_dir,'evaluation_matrix','alpha_',num2str(alpha),'beta_',num2str(beta),'\'];
evaluation_matrix_savepath = [evaluation_matrix_savedir,type,'_',name,'.mat'];
evaluation_matrix = stats;
if exist(evaluation_matrix_savedir,'dir')==0
    mkdir(evaluation_matrix_savedir)
end
 save(evaluation_matrix_savepath,'evaluation_matrix')   


%  total_time =  toc(t1);
%  display(['Timings for brain: ',type,'_',name,':']) 
%  display(['Total processing: ',num2str(total_time) ])
%  display(['kkn: ',num2str(knn_time)]);
%  display(['graphcut: ',num2str(graphcut_time)])
%  total_score = total_score + avg;
% save('all_variables.mat');
% clear all
% load('all_variables.mat');

 %   end
 end
 % record_method_evaluation(alpha, beta, save_path3, main_result_dir);
%  parameter_vector = [alpha, beta, total_score/30]';
%  parameters = [parameters parameter_vector'];

    end
end
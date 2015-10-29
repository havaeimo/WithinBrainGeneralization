clear all
close all;

mex -lut /home/local/USHERBROOKE/havm2701/git.repos/semi_bts_knn/mex_test/gcmex-2.3.0/GCMex/markov_network_multilable_general3Dgraph.cpp /home/local/USHERBROOKE/havm2701/git.repos/semi_bts_knn/mex_test/gcmex-2.3.0/GCMex/GCoptimization.cpp /home/local/USHERBROOKE/havm2701/git.repos/semi_bts_knn/mex_test/gcmex-2.3.0/GCMex/graph.cpp /home/local/USHERBROOKE/havm2701/git.repos/semi_bts_knn/mex_test/gcmex-2.3.0/GCMex/LinkedBlockList.cpp /home/local/USHERBROOKE/havm2701/git.repos/semi_bts_knn/mex_test/gcmex-2.3.0/GCMex/maxflow.cpp
 mex -lut /home/local/USHERBROOKE/havm2701/git.repos/semi_bts_knn/mex_test/gcmex-2.3.0/GCMex/interactiveGraphcut_withdataterms.cpp /home/local/USHERBROOKE/havm2701/git.repos/semi_bts_knn/mex_test/gcmex-2.3.0/GCMex/GCoptimization.cpp /home/local/USHERBROOKE/havm2701/git.repos/semi_bts_knn/mex_test/gcmex-2.3.0/GCMex/graph.cpp /home/local/USHERBROOKE/havm2701/git.repos/semi_bts_knn/mex_test/gcmex-2.3.0/GCMex/LinkedBlockList.cpp /home/local/USHERBROOKE/havm2701/git.repos/semi_bts_knn/mex_test/gcmex-2.3.0/GCMex/maxflow.cpp
global T1C T2 T1 FLAIR
classifier_name = 'kNN';
user_idx = 'test2red_EDIT_3';
spatial_feature = 1; % set to one if you want the spatial features
k=3;
compute_statistics =0;
if ~strcmp(classifier_name,'kNN')
    k = NaN;
end
params = [0.033,0.004];

param = params;
alpha = param(1);
beta = param(2);



% output directories
if spatial_feature==0
    result_sub_directory = '3dimensional_F_newmasks';
else
    result_sub_directory = '6dimensional_F_newmasks';
end




if compute_statistics ==0
    main_result_dir = ['/home/local/USHERBROOKE/havm2701/data/RESULTS/within_brain_classification/BRATS_2015_redbrains_test2_res',filesep,user_idx,filesep,'challenge_results',filesep,'',result_sub_directory,'',filesep,'',classifier_name,'_alpha_',num2str(alpha),'_beta_',num2str(beta),'',filesep,''];
elseif compute_statistics==1
    main_result_dir = ['/home/local/USHERBROOKE/havm2701/data/RESULTS/within_brain_classification/BRATS_2015_CHALLENGE',filesep,'journal_2014_linearSVM',filesep,'training_results',filesep,'',result_sub_directory,'',filesep,'',classifier_name,'_alpha_',num2str(alpha),'_beta_',num2str(beta),'',filesep,''];
end

%% define paths and parameters
Nb_Selected = 2000;

Nb_healthy = 2000;
downsamplerate = 20;


knn_result_save_dir =[main_result_dir,classifier_name,'',filesep,''];
if exist(knn_result_save_dir,'dir')==0
    mkdir(knn_result_save_dir)
end


knn_medianfilter_result_save_dir =[main_result_dir,classifier_name,'_medianfilter',filesep,''];
if exist(knn_medianfilter_result_save_dir,'dir')==0
    mkdir(knn_medianfilter_result_save_dir)
end

boykov_result_save_dir =  [main_result_dir,'boykov',filesep,'nodenoising',filesep,''];
if exist(boykov_result_save_dir,'dir')==0
    mkdir(boykov_result_save_dir)
end





boykov_denoised_result_save_dir =  [main_result_dir,'boykov',filesep,'Denoised',filesep,''];
if exist(boykov_denoised_result_save_dir,'dir')==0
    mkdir(boykov_denoised_result_save_dir)
end


    brain_dir = ['/home/local/USHERBROOKE/havm2701/data/Data/BRATS/BRATS2015/brats15_test/'];
    mask_dir =  ['/home/local/USHERBROOKE/havm2701/data/RESULTS/within_brain_classification/redtest1_porcessedMasks/theproblematicbrain'];

   
mask_list = dir(mask_dir);
 mask_list = mask_list(3:end);
mask_type_path = mask_dir;
% mask_type_list = mask_type_list(3:end);
%  for subfoldercounter=1:length(mask_type_list)
%      type = mask_type_list(subfoldercounter).name;
%      mask_type_path = [mask_dir,type];
%      mask_list = dir(mask_type_path);
%      mask_list = mask_list(3:end);
type = 'HL';
     for brain_counter =1:length(mask_list)
          
         name = mask_list(brain_counter).name;
         mask_name_path = [mask_type_path,'',filesep,'',name];
         brain_name_path = [brain_dir,type,'_',name];

t1 = tic();
%% load modalities
compute_statistics=0;
        [T1C, T2, FLAIR, truth,T1,info,flair_name]  = load_modalities(brain_name_path,compute_statistics);
        [height width depth] = size(T1C);
        T1C = double(T1C);
        T2 = double(T2);
        FLAIR = double(FLAIR);
        if compute_statistics == 0
            truth = zeros(size(T1C));
        end
        % take backup of the modalities
        truth_backup = truth;

%% reduce the resolution of T2 (check for Maxime)
%        T2 = rescalingT2(T2);
        
%% creating the mask and ROI
   % MASK = make_mask_b15(T1C, T2, T1);
% or load the mask from a path
%mask_name_path = '/media/Expansion Drive/User/Research/Results/BRATS-2_KNN_ENDAUGUST',filesep,'HG',filesep,'0026';
 %mask_name_path = ['/media/Expansion Drive/User/Research/Results/BRATS-2_KNN_ENDAUGUST',filesep,'',type,'',filesep,'',name];
 %mask_name_path = [mask_dir,type,'',filesep,'',name];
 MASK = load_mask(mask_name_path);
 backup_MASK = MASK;
% you can Edit the mask
     %MASK = edit_mask_b15(T1C, T2, T1, MASK);
 
%[ymin,ymax,xmin,xmax,zmin,zmax] = make_regtangle(T1C);  
% save ROI (regtangle) coordinate points
 


%% Detect the borders of the brain

[ymin,ymax,xmin,xmax,zmin,zmax] = find_boundingbox_borders(T1C);


 
disp('***********************************')
disp(['Processing brain ',type,'_',name, ' with k = ', num2str(k)])
disp('***********************************')


%%
T1C = T1C(xmin:xmax , ymin:ymax , zmin:zmax);
T2 = T2(xmin:xmax , ymin:ymax , zmin:zmax);
T1 = T1(xmin:xmax , ymin:ymax , zmin:zmax);
FLAIR = FLAIR(xmin:xmax , ymin:ymax , zmin:zmax);

MASK = MASK(xmin:xmax , ymin:ymax , zmin:zmax);

%%  make space and selected points space( MADE CHANGES FOR MAKING TEXT FILES >> FIX LATER)
[space,mask_idx,truth_idx] = make_space(T1C, FLAIR, T2,truth, MASK);

 space2 = space;
         background = find(sum(space(1:3,:))< 0.001); 
         space(:,background) = [];
         mask_idx(background) = [];
         if compute_statistics
         truth_idx(background) = [];
         end




[selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );


  t2 = tic() ;
 [h,w,d] = size(T1C);
%% classifiers

 
% KNN search classifier
   
    kkk=k;

    [IDX,D] = knnsearch(selected_space',space','K', k ,'NSMethod','kdtree' ,'Distance','euclidean');

    % make data-term_KNN
    num_classes = 4;

    matrix = make_dataterm_matrix(IDX,mask_idx,space,T1C,num_classes);
    matrix = (matrix/k);

    energy_base = 1 - matrix;



%% APPLY graphcut BOYKOV JOLLY with dataterms
% alpha = 0.001;


 energy = alpha * energy_base;
%  beta = .0053;
t3 = tic();
 lf = interactiveGraphcut_withdataterms(space2,double(MASK),energy,beta);
time_graphcut = toc(t3);
 [h,w,d] = size(T1C);
  lf = reshape(lf,h, w, d);
   segmented_volume = lf ;

% segmented_brain = zeros(height,width,depth);
% segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
% 
%  segmented_brain(segmented_brain == 3)=4;
%  segmented_brain(segmented_brain == 1)=3;
% 
% %%%%%%%%%%%% save results
% mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);
% mha_result_path = [boykov_result_save_dir,mha_result_name];
% segmented_brain = uint16(segmented_brain);
% writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset); 

   
t4 = tic();
% apply median filter
segmented_volume = medfilt3(segmented_volume, [3,3,3]);
% compute statistics 
time_medianfilter = toc(t4);
segmented_brain = zeros(height,width,depth);
segmented_brain(xmin:xmax , ymin:ymax , zmin:zmax) = segmented_volume;
truth = truth_backup;
 segmented_brain(segmented_brain == 3)=4;
 segmented_brain(segmented_brain == 1)=3;
 segmented_brain(segmented_brain == 3)=1;
%%%%%%%% save results 
mha_result_name = strrep(flair_name , 'Brain.XX.O.MR_Flair',['Seg_',type,'_',name]);

mha_result_path = [boykov_denoised_result_save_dir,mha_result_name];
segmented_brain = uint16(segmented_brain);
writemetaimagefile(mha_result_path, segmented_brain, info.PixelDimensions,info.Offset); 
     end
%end
 


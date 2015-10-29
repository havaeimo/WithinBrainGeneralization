clear all
close all
global T1C T2 T1 FLAIR

%% define paths and parameters
Nb_Selected = 500;
k=21;
Nb_healthy = 500;
downsamplerate = 20;
%results_path_root ='/home/local/USHERBROOKE/havm2701/data/RESULTS/within_brain_classification/BRATS_2015_Testpip/';
%new_results_path_root ='/home/local/USHERBROOKE/havm2701/data/RESULTS/within_brain_classification/BRATS_2013_LEADERBOARD_MASK2/';
% results_path_root =['/home/local/USHERBROOKE/havm2701/data/RESULTS/within_brain_classification/BRATS_2015_challengeday','_2classMASKS'];
% if exist(results_path_root,'dir')==0
%     mkdir(results_path_root)
% end
% main_result_dir = results_path_root;

%results_path_root = '/Users/uoft/Dropbox/PhD/BratsAnalysis/CHALLENGE2013_KNN_mha/';
% type = '';
% name='brats_tcia_pat383_0001';
%  name_hl = ['HL_',name];
% % Brain_path = ['/home/local/USHERBROOKE/havm2701/data/Data/BRATS/BRATS2015/brats15_test2','/',name_hl];
% 
%  Brain_path = ['/home/local/USHERBROOKE/havm2701/data/Data/BRATS/BRATS2015/brats15_test_preprocess/BRATS15_TEST2','/',name_hl];
% %Brain_path = ['/home/local/USHERBROOKE/havm2701/data/Data/BRATS/BRATS2015/BRATS_BRAINS_2015_preprocessed/BRATS2014_training','/',type,'/',name_hl];
% % Brain_path = ['/Volumes/FreeAgent Drive/Research UdeSh/BRATS-2/Image_Data/',type,'/',name];
%  %%%mask_path = ['N:\User\Research\Results\BRATS-2_KNN_ENDAUGUST\',type,'\',name]
% mask_path = [results_path_root,'/',name] 
%  if exist(results_path_root,'dir')==0
%     mkdir(results_path_root)
% end
% result_path_brain = [results_path_root,'/', name];
% if exist(result_path_brain,'dir')==0
%     mkdir(result_path_brain)
% end
 type = 'HL';
 name='brats_2013_pat0135_1';
 name_hl = ['HL_',name];
 %Brain_path = ['/home/local/USHERBROOKE/havm2701/data/Data/BRATS/BRATS2015/brats15_test','/',name_hl];
Brain_path = ['/home/local/USHERBROOKE/havm2701/data/Data/BRATS/BRATS2015/brats15_test2','/',name_hl];
%Brain_path = ['/home/local/USHERBROOKE/havm2701/data/Data/BRATS/BRATS2015/brats15_test_preprocess/BRATS15_TEST2','/',name_hl];

 % Brain_path = ['/Volumes/FreeAgent Drive/Research UdeSh/BRATS-2/Image_Data/',type,'/',name];
mask_path = ['/home/local/USHERBROOKE/havm2701/data/RESULTS/within_brain_classification/BRATS_2015_challengeday_2class_2013_MASKS/',name];


%% load modalities
compute_statistics = 0;
[T1C, T2, FLAIR, truth,T1,info,flair_name] = load_modalities(Brain_path,compute_statistics);
[height width depth] = size(T1C);
T1C = double(T1C);
T2 = double(T2);
T1 = double(T1);

%% creating the mask
%MASK = make_mask(T1C, T2, T1);
% or load the mask from a path (future work)
%mask_path = result_path_brain;
%MASK = load_mask(mask_path);
MASK = zeros(size(T1C));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%MAKE MASK%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 display('***** Make New Mask*****')
[height,width,depth] = size(T1C);

while 1
    %MASK_slice = imread('Mask_HG14_slice80.png');
    while 1
        fig2=figure();
         d=input('press "1" to select healthy region');
        if d==1
            sliceH = input('enter slice number for healthy (from image J):');
%             sliceH = depth - sliceH;
            %modality = input('Enter modality to slelect from: ','s');
            modality = 'FLAIR';
            view = input('Enter view to slelect from: ');
            imshow_modality(sliceH,modality,view);
            
            roi= roipoly;
            if view ~= 3 
                roi = flipud(roi)';
            elseif view ==3
                roi = roi';
            end
            MASK_sliceH = roi;
            MASK_sliceH = 10 * MASK_sliceH;
            imshow(roi);
            roi_index = find(roi>0);

            if view ==1
            MASK_sliceH =MASK(sliceH,:,:);
            MASK_sliceH(roi_index)= 10;
            MASK(sliceH,:,:) = MASK_sliceH;
            elseif view ==2
            MASK_sliceH =MASK(:,sliceH,:);
            MASK_sliceH(roi_index)= 10;
            MASK(:,sliceH,:) = MASK_sliceH;
            elseif view == 3
            MASK_sliceH =MASK(:,:,sliceH);
            MASK_sliceH(roi_index)= 10;
            MASK(:,:,sliceH) = MASK_sliceH;
            end
        else 
            break
        end                
    end
    while 1
        d=input('press "1" to select edema region');
        if d==1
            sliceE = input('enter slice number for edema (from image J):');
%             sliceE = depth - sliceE;
            %modality = input('Enter modality to slelect from: ','s');
           modality = 'FLAIR';
            view = input('Enter view to slelect from: ');
            imshow_modality(sliceE,modality,view);
             roi= roipoly;
            if view ~= 3 
                roi = flipud(roi)';
            elseif view ==3
                roi = roi';
            end
            imshow(roi);
            roi_index = find(roi>0);

            if view ==1
            MASK_sliceE =MASK(sliceE,:,:);
            MASK_sliceE(roi_index)= 2;
            MASK(sliceE,:,:) = MASK_sliceE;
            elseif view ==2
            MASK_sliceE =MASK(:,sliceE,:);
            MASK_sliceE(roi_index)= 2;
            MASK(:,sliceE,:) = MASK_sliceE;
            elseif view == 3
            MASK_sliceE =MASK(:,:,sliceE);
            MASK_sliceE(roi_index)= 2;
            MASK(:,:,sliceE) = MASK_sliceE;
            end
        else 
            break
        end                
    end

   

    dd=input('press "1" to continue selecting ROI');
     if dd~=1
          break
     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 MASK = double(MASK);
%MASK = edit_mask(T1C, T2, T1,MASK);





% %% assign labels of the mask form truth
% % putting the healthy labels to 10
% h_index = find(truth==0);
% truth(h_index) = 10;
% 
% 
%  Background_index = find(MASK > 0);
%  MASK(Background_index) = truth(Background_index);
 nicroses_idx = find(MASK ==1);
 MASK(nicroses_idx) = 3; 

%% make space and selected points space
% [space,mask_idx] = make_space(T1C, T1, T2,truth, MASK);
% 
% [selected_space, mask_idx] = make_selected_space(space, mask_idx, Nb_healthy,Nb_Selected );

%% KNN search

% importance_map = KNN_search(space, selected_space, k);
% [IDX,D] = KNN_search(selected_space,space, k)
% [IDX,D] = knnsearch(selected_space',space','K', k ,'NSMethod','kdtree' ,'Distance','euclidean');

%% assign labels
%segmented = assign_labels(IDX,D,mask_idx, space, height, width, depth);

%t=tic;

%% apply median filter
%segmented = medfilt3(segmented, [5,5,5]);
%   segmented = u_medfilt3(segmented,5);
%time = toc(t);
%% compute statistics 
% [dice,jaccard,precision,recall] = compute_eval_multilabel_metrics(MASK,truth);


%% save results
% save the mha file for the mask
%MASK = uint16(MASK);
%mha_MASK_path = [mask_path,'/','MASK_',type,'_',name,'.mha'];
%writemetaimagefile(mha_MASK_path, MASK, info.PixelDimensions,info.Offset);
MASK = uint16(MASK);
  if exist(mask_path,'dir')==0
     mkdir(mask_path)
 end
mha_MASK_path = [mask_path,'/','MASK_','_',name,'.mha'];
writemetaimagefile(mha_MASK_path, MASK, info.PixelDimensions,info.Offset);

%save nii files
%mask_nii = make_nii(MASK);
%save_nii(nii, [result_path_brain,'\','6dimportancemap.nii.gz']);
%save_nii(mask_nii, [result_path_brain,'\','MASK.nii.gz']);
%segmented_nii = make_nii(segmented); 
%save_nii(segmented_nii, [result_path_brain,'\',type,'_',name,'_segmented_tumor.nii.gz']);

%%save statistics
%save([result_path_brain,'\',type,'_',name,'_stat_results.mat'],'precision','recall');

%%save meta files(.mha)
%mha_result_name = strrep(T1_name , 'XX.O.MR_T1',['Seg_',type,'_',name]);
%mha_result_path = [result_path_brain,'\',mha_result_name];
%segmented=uint16(segmented);
%writemetaimagefile(mha_result_path, segmented, info.PixelDimensions,info.Offset);

clear all
close all
Nb_Selected = 500;
downsamplerate = 20;
name = 'BRATS_LG0004';
dir_name =['J:\User\Research\Brain_data\BRATS-1\Images\',name];

results_path_root = 'C:\Users\havm2701\Dropbox\PhD\BratsAnalysis\CHALLENGE2013\';
if exist(results_path_root,'dir')==0
    mkdir(results_path_root)
end
result_path_brain = [results_path_root, name];
if exist(result_path_brain,'dir')==0
    mkdir(result_path_brain)
end
    %dir_name =['C:\Users\havm2701\Documents\MATLAB\',name];

list_modality = dir (dir_name);
list_modality = list_modality(2:end);
for j=1:length(list_modality)
    name_modality = list_modality(j);
    if ~isempty(strfind(name_modality.name,'T1C.nii'))
        filepath = [dir_name,'\',name_modality.name];
        T1C = load_nii(filepath, [], 1);
        T1C=T1C.img; 
    end
    if ~isempty(strfind(name_modality.name,'FLAIR.nii'))
        filepath = [dir_name,'\',name_modality.name];
        FLAIR = load_nii(filepath, [], 1);
        FLAIR=FLAIR.img; 
    end
    if ~isempty(strfind(name_modality.name,'T2.nii'))
        filepath = [dir_name,'\',name_modality.name];
        T2 = load_nii(filepath, [], 1);
        T2=T2.img; 
    end
    if ~isempty(strfind(name_modality.name,'truth.nii'))
        filepath = [dir_name,'\',name_modality.name];
        truth = load_nii(filepath, [], 1);
        truth=truth.img; 
    end
end

[height width depth] = size(T1C);

%     T1C= T1C(:);
%     T2 = T2(:);
T1C = double(T1C);
T2 = double(T2);
FLAIR = double(FLAIR);
T1C = T1C/max(max(max(T1C)));
T2 = T2/max(max(max(T2)));
FLAIR = FLAIR/max(max(max(FLAIR)));
truth(truth<0.5)=0;                 %healthy 
truth(0.5<truth & truth<1.5)=2;     %tomur
truth(1.5<truth & truth<2.5)=1;     %edema


%% creating the mask


MASK = zeros(size(T1C));
%MASK_slice = imread('Mask_HG14_slice80.png');
decide = false;
while decide == false
    fig2=figure();
    d=input('press "1" to select edema region');
    if d==1
        sliceE = input('enter slice number for edema (from image J):');
        sliceE = depth - sliceE;
        imshow(T1C(:,:,sliceE),[])
        MASK_sliceE = roipoly;
        MASK_sliceE = 2 * MASK_sliceE;
        imshow(MASK_sliceE);
        MASK(:,:, sliceE)= min((MASK(:,:, sliceE) +MASK_sliceE),2);

    end
     d=input('press "1" to select enhancing tumor');
    if d==1
        sliceET = input('enter slice number for enhancing tumor:');
        sliceET = depth - sliceET;
        imshow(T1C(:,:,sliceET),[])
        MASK_sliceET = roipoly;
        MASK_sliceET = 4* MASK_sliceET;
        imshow(MASK_sliceET);
        MASK(:,:, sliceET)= min((MASK(:,:, sliceET) +MASK_sliceET),4);
    end
    d=input('press "1" to select none_enhancing tumor');
    if d==1
        sliceNT = input('enter slice number for non enhancing tumor:');
        sliceNT = depth - sliceNT;
        imshow(T1C(:,:,sliceNT),[])
        MASK_sliceNT = roipoly;
        MASK_sliceNT = 3* MASK_sliceNT;
        imshow(MASK_sliceNT);
        MASK(:,:, sliceNT)= min((MASK(:,:, sliceNT) +MASK_sliceNT),3);

    end
    
    d=input('press "1" to continue selecting ROI');
    if d~=1
        break
    end
end

%     imshow
%     MASK_slice(MASK_slice>0)=1;

%% creating the 6 D space (T1,T2,FLAIR,x,y,z)
count = 1;
space = zeros(6,length(T1C(:)));
mask_idx = zeros(1,size(space,2));
% [width, height, depth] = size(T1C);
 for DEP = 1:depth
    for ROW = 1:height
        for COL=1:width  
        space(:,count)=[T1C(ROW,COL,DEP);T2(ROW,COL,DEP);FLAIR(ROW,COL,DEP);ROW/height;COL/width;DEP/depth];
        mask_idx(count) = MASK(ROW,COL,DEP);
        count = count+1 ;
        end
    end
end

%     for i=1:size(space,1)
%         space(i,:) = space(i,:)/max(space(i,:));
%     end
% idx_downsample = downsample(1:size(space,2), downsamplerate);
% space = space(:,idx_downsample);
selected_points = find(mask_idx>0); 
iindex = downsample(selected_points, floor(length(selected_points)/Nb_Selected)+1);%sort(uint32(rand(Nb_Selected,1)*length(selected_points)));
selected_space = space(:,iindex);
mask_idx = mask_idx(iindex);

edema_points = find(mask_idx==2);
E_tumor_points = find(mask_idx==4);
N_tumor_points = find(mask_idx==3);
%iindex = sort(uint32(rand(Nb_Selected,1)*length(selected_points))); %select some random points from the selected ROI
%iindex = myunique(iindex);  %remove repeated indexes
selected_spaceET = selected_space(:,E_tumor_points);
selected_spaceNT = selected_space(:,N_tumor_points);

selected_spaceE = selected_space(:,edema_points);

% space = space(:,:);   
% space_repmat = repmat(space,[1 1 size(selected_space,2)]);
% selected_space_3d = reshape(selected_space,size(selected_space,1),1,size(selected_space,2));
% selected_space_repmat = repmat(selected_space_3d,[1 size(space,2) 1]);
importance_map = zeros (2,size(space,2));
class = [];
for i = 1:size(space,2)
    
point = space(:,i);
if ~isempty(selected_spaceET);
    tmpET = abs(selected_spaceET-repmat(point,1,size(selected_spaceET,2)));
    tmpET = sum(tmpET,1);
    valET = min(tmpET);
    class = [class , [valET; 4]];
end
if ~isempty(selected_spaceNT)
    tmpNT = abs(selected_spaceNT-repmat(point,1,size(selected_spaceNT,2)));
    tmpNT = sum(tmpNT,1);
    valNT = min(tmpNT);
    class = [class,[valNT;3]];
end
if ~isempty(selected_spaceE)
    tmpE = abs(selected_spaceE-repmat(point,1,size(selected_spaceE,2)));
    tmpE = sum(tmpE,1);
    valE = min(tmpE);
    class = [class,[valE;2]];
end
[val,index_class] = min(class(1,:));

importance_map(1,i) = val;
importance_map(2,i) = class(2,index_class) ;      %signifiying the voxel is more likely to be edema


end

%        for j = 1:size(selected_space,2)
%            p_selected = selected_space(1:5,j);
%            distance = sqrt(sum((p_selected - point).^2));
%            imp = min(imp,distance);
%        end
%        importance_map(i) = imp;
%    end

importance_nii = zeros(size(T1C));   
importance_nii_idx = zeros(size(T1C));   

count = 1;
for DEP = 1:depth
    for ROW = 1:height
        for COL=1:width

        importance_nii(ROW,COL,DEP) = importance_map(1,count);
        importance_nii_idx(ROW,COL,DEP) = importance_map(2,count);  
        count = count+1 ;
        end
    end
end   
nii = make_nii(importance_nii);
mask_nii = make_nii(MASK);
save_nii(nii, [result_path_brain,'\','6dimportancemap.nii.gz']);
save_nii(mask_nii, [result_path_brain,'\','MASK.nii.gz']);
 


fig5=figure(5);
pr = -Inf;
stat_result=[];
for t = 0.1:0.005:0.5
segmented_nii = importance_nii;    
segmented_nii(segmented_nii > t) =1;
segmented_nii(segmented_nii<=t) = 0;


segmented_nii = ceil(-(segmented_nii - 0.5));
idx = find(segmented_nii>0);
segmented_nii(idx) = importance_nii_idx(idx);
[dice,jaccard,precision,recall] = compute_eval_metrics(segmented_nii,truth);
plot(precision,recall,'.b')
xlabel('precision')
ylabel('recall')
title(['precision_recall for threshold change on brain ',name])
axis([0 100 0 100])
if precision+recall > pr
pr = precision+recall;
topt = t;
precision_opt = precision;
recall_opt = recall;
end
hold on
result=[t;precision;recall];
stat_result = [stat_result,result];
end

plot(precision_opt,recall_opt,'.green')
hold off
saveas(fig5, [result_path_brain,'\','6Dmetod_PrecisionRecall.png']);
segmented_nii = importance_nii;    
segmented_nii(segmented_nii > topt) =1;
segmented_nii(segmented_nii<=topt) = 0;
segmented_nii = ceil(-(segmented_nii - 0.5));
segmented_nii = medfilt3(segmented_nii,5);
idx = find(segmented_nii>0);
segmented_nii(idx) = importance_nii_idx(idx);
[dice,jaccard,precision,recall] = compute_eval_metrics(segmented_nii,truth);
segmented_nii = make_nii(segmented_nii);    
save_nii(segmented_nii, [result_path_brain,'\',name,'_segmented_tumor.nii.gz']);
save([result_path_brain,'\',name,'_stat_results.mat'],'topt','stat_result','precision','recall');
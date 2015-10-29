function [space,mask_idx,truth_idx] = make_space(T1C, FLAIR, T2,truth, MASK) 



T1C = double(T1C);
T2 = double(T2);
FLAIR = double(FLAIR);
T1C = T1C/max(max(max(T1C)));
T2 = T2/max(max(max(T2)));
FLAIR = FLAIR/max(max(max(FLAIR)));


% 
% truth(truth<0.5)=0;                 %healthy 
% truth(0.5<truth & truth<1.5)=2;     %tomur
% truth(1.5<truth & truth<2.5)=1;     %edema

%% creating the 6 D space (T1,T2,FLAIR,x,y,z)

space     = zeros(6,length(T1C(:)));
mask_idx  = zeros(1,size(space,2));
truth_idx = zeros(1,size(space,2));
[height,width, depth] = size(T1C);
i=1; 
% IMPORTANT : it is essential that the matrix space be made in F-order(column wise).
% this is how matlab stors the matrices for mex files. 
 for DEP = 1:depth
    for COL=1:width
        for ROW = 1:height
            
                space(:,i)=[T1C(ROW,COL,DEP);T2(ROW,COL,DEP);FLAIR(ROW,COL,DEP);ROW/height;COL/width;DEP/depth];
                mask_idx(i) = MASK(ROW,COL,DEP);
               
                i = i+1;
                
        end
    end
end


%make text file moved to the main m file
if isnan(truth)
    truth_idx = nan;
else 
  for DEP = 1:depth
    for COL=1:width
        for ROW = 1:height
                truth_idx(i) = truth(ROW,COL,DEP);
                i = i+1;
                
        end
    end
  end
end

%mask_idx(background) = 10;
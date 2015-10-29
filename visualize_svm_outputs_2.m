function [segmented posterior_matrix] = visualize_svm_outputs_2(inputfile,info,probability_flag,backup_MASK,nb_classes,outputfile)

% this method converts a text file to a mha file. the resulting mha file
% can be visualised with ITK-snap.
% the hight width depth dimensions of the 3d matrix must be given to the
% method
%input arguments:
%   inputfile : path to the input file
%   outputfile: path to the output file
%   info: header information for the mha file
%   probability_flag: =1 if the input textfile contains probability information    
% october 2013
if nargin < 5
   outputfile = nan;
   nb_classes = 4;
   produce_segmented_image = 0; 
   
end
if ~isnan(outputfile)
    produce_segmented_image = 1;
end

classes = unique (backup_MASK(:));
classes([1,end],:)= classes([end,1],:);
classes(end) = [];
nb_classes = length(classes); 
size = info.Dimensions ;
height = size(1);
width  = size(2);
depth  = size(3);
posterior_matrix = zeros(nb_classes,height,width,depth);
%Posterior_matrix(1,:,:,:)=1;
segmented = zeros(height, width, depth);
f = fopen(inputfile,'r');
line = fgetl(f);
while ischar(line)
    if  probability_flag == 0
        [label x y z] = read_libsvm_line(line);
        ROW = uint8(x * height);
        COL = uint8(y * width);
        DEP = uint8(z * depth);
        segmented(ROW,COL,DEP) = label;
        line = fgetl(f);
    elseif probability_flag == 1
         [label,spatial_info, a] = read_libsvm_probability_line_2(line,classes);
        ROW = uint8(spatial_info(1) * height);
        COL = uint8(spatial_info(2) * width);
        DEP = uint8(spatial_info(3) * depth);

        line = fgetl(f);

        for id = 1: nb_classes
            posterior_matrix(id,ROW,COL,DEP) = a(id);
        end
        
    end

end



for c=1:length(classes)
    idx_mask_c = find(backup_MASK == classes(c));
    ix = randperm(length(idx_mask_c));
    random_point_c = idx_mask_c(ix(1));
    [dummy,c_id] = max(posterior_matrix(:,random_point_c));
    posterior_matrix([c,c_id],:,:) = posterior_matrix([c_id,c],:,:);
end

zero_idx = find(sum(posterior_matrix,1)==0);
posterior_matrix(1,zero_idx) = 1;

for i = 1:height
    for j = 1: width
        for k = 1: depth
            [dummy, label] = max(posterior_matrix(:,i,j,k));
            segmented(i,j,k)= label;
        end
    end
end
        

fclose(f);
segmented(segmented==1)=0;

if produce_segmented_image == 1
    segmented=uint16(segmented);
    writemetaimagefile(outputfile, segmented, info.PixelDimensions,info.Offset);
end


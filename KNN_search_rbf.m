function [IDX,D] = KNN_search_rbf(selected_space, space, k,sig)

% space = space(:,:);   
% space_repmat = repmat(space,[1 1 size(selected_space,2)]);
% selected_space_3d = reshape(selected_space,size(selected_space,1),1,size(selected_space,2));
% selected_space_repmat = repmat(selected_space_3d,[1 size(space,2) 1]);
selected_space = selected_space';
space = space';

IDX = zeros (size(space,2),k);
D = zeros (size(space,2),k);
for i = 1:size(space,2)
    
point = space(:,i);
tmp = ((selected_space-repmat(point,1,size(selected_space,2))).^2);
tmp1 = tmp(1:3,:);
tmp2 = tmp(4:6,:);
d = 1- gaussmf(tmp2,[sig 0]) ;
tmp = [tmp1;d];


tmp = sum(tmp,1);
[tmp IX] = sort(tmp);
D(i,:) = tmp(1:k);
IDX (i,:)= IX(1:k);



end

function [IDX,D] = KNN_search_polinomial(selected_space, space, kkk,m)

if size(space,1)<size(space,2)
    error('Error :dimension mis-match')
end

projected_selected_space= [];
projected_space= [];
for i = 1: m
    temp_selected_space = selected_space.^m;
    projected_selected_space = [projected_selected_space temp_selected_space];
    
    temp_space = space.^m;
    projected_space = [projected_space temp_space];
end
    
[IDX,D] = knnsearch(projected_selected_space,projected_space,'K', kkk  ,'Distance','euclidean');    
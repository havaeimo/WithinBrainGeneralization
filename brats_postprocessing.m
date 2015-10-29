clear all
close all


%% define paths and parameters

Ref_path = '/home/local/USHERBROOKE/havm2701/Desktop/Matlab_brats/ref'
Query_path = '/home/local/USHERBROOKE/havm2701/Desktop/Matlab_brats/query'

ref_list = dir(Ref_path)
query_list = dir(Query_path)

ref_list = ref_list(3:end)
query_list = query_list(3:end)

for i = 1:length(query_list)
   if length(strfind(ref_list(i).name,'p')) ==0
       continue
   end
   id = strfind(query_list(i).name,'_')
   id  = id(1);
   query_name = query_list(i).name(id+1:id+4)
   display query_name

   for j = 1:length(ref_list)
   if isempty(strfind(ref_list(j).name,'process'))
       continue
   end
   
   if isempty(strfind(ref_list(j).name,query_name ))
       continue
   end
  
   query_seg_path = [Query_path,'/',query_list(i).name];
   ref_seg_path = [Ref_path,'/',ref_list(j).name];
   query_seg = uint16(mha_read_volume(query_seg_path));
   info  = mha_read_header(query_seg_path);
   ref_seg = uint16(mha_read_volume(ref_seg_path));
   idx = (ref_seg > 0);
   idy = (query_seg > 0);
   mask = uint16(xor(idx,idy));
   
   idx = uint16(idx);  
   new_query_seg = query_seg .* idx;
   new_query_seg = (new_query_seg);
   result_name = ['/home/local/USHERBROOKE/havm2701/Desktop/Matlab_brats/Results/',query_list(j).name];
   result_name_mask = ['/home/local/USHERBROOKE/havm2701/Desktop/Matlab_brats/Masks/Mask_',query_list(j).name];

   writemetaimagefile(result_name, new_query_seg, info.PixelDimensions,info.Offset);
   writemetaimagefile(result_name_mask, mask, info.PixelDimensions,info.Offset);
   
   
   end
end    

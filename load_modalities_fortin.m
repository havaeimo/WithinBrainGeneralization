function [T1C, T2 , T1] = load_modalities_fortin(dir_name)


list_modality = dir (dir_name);
list_modality = list_modality(2:end);
for j=1:length(list_modality)
    name_modality = list_modality(j);
    if ~isempty(strfind(name_modality.name,'t1_gad'))
        filepath = [dir_name,'\',name_modality.name];
        T1C = load_nii(filepath, [], 1);
        T1C=T1C.img; 
    end

    if ~isempty(strfind(name_modality.name,'b0'))
        filepath = [dir_name,'\',name_modality.name];
        T2 = load_nii(filepath, [], 1);
        T2=T2.img; 
    end
    
        if ~isempty(strfind(name_modality.name,'t1_final'))
        filepath = [dir_name,'\',name_modality.name];
        T1 = load_nii(filepath, [], 1);
        T1=T1.img; 
    end

end

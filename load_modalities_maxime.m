function [T1C, T2, T1,T1_C_info] = load_modalities_maxime(Brain_path)
%% load 4 modalities
file_list = dir(Brain_path);
file_list = file_list(3:end);

for k = 1:length(file_list)
             filename = file_list(k);
             if ~isempty(strfind(filename.name,'T1')) && isempty(strfind(filename.name,'T1c'))
                filepath = [Brain_path,filesep,filename.name];
                files = dir(filepath);
                files = files(3:end);
                for kk = 1: length(files)
                    modality=files(kk);
                    if ~isempty(strfind(modality.name,'T1'))
                        modality_path =  [Brain_path,filesep,filename.name,filesep,modality.name];
                        T1 = mha_read_volume(modality_path);
                        info  = mha_read_header(modality_path);
                    end
                end 
            end
  
            if ~isempty(strfind(filename.name,'T2'))
                filepath = [Brain_path,filesep,filename.name];
                files = dir(filepath);
                files = files(3:end);
                for kk = 1: length(files)
                    modality=files(kk);
                    if ~isempty(strfind(modality.name,'T2'))
                        modality_path =  [Brain_path,filesep,filename.name,filesep,modality.name];
                        T2 = mha_read_volume(modality_path);
                        info  = mha_read_header(modality_path);
                    end
                end 
            end  
%             if ~isempty(strfind(filename.name,'Flair'))
%                 filepath = [Brain_path,filesep,filename.name];
%                 files = dir(filepath);
%                 files = files(3:end);
%                 for kk = 1: length(files)
%                     modality=files(kk);
%                     if ~isempty(strfind(modality.name,'Flair'))
%                         flair_name = modality.name;
%                         modality_path =  [Brain_path,filesep,filename.name,filesep,modality.name];
%                         FLAIR = mha_read_volume(modality_path);
%                         info  = mha_read_header(modality_path);
%                     end
%                 end 
%             end  
            if ~isempty(strfind(filename.name,'T1c'))
                filepath = [Brain_path,filesep,filename.name];
                files = dir(filepath);
                files = files(3:end);
                for kk = 1: length(files)
                    modality=files(kk);
                    if ~isempty(strfind(modality.name,'T1c'))
                        modality_path =  [Brain_path,filesep,filename.name,filesep,modality.name];
                        T1C = mha_read_volume(modality_path);
                        T1_C_info  = mha_read_header(modality_path);
                    end
                end 
            end  

end      
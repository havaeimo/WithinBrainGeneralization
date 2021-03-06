function [T1C, T2, FLAIR, truth,T1,info_flair,flair_name] = load_modalities(Brain_path,compute_statistics)
%% load 4 modalities
file_list = dir(Brain_path);
file_list = file_list(3:end);

for k = 1:length(file_list)
             filename = file_list(k);
            if ~isempty(strfind(filename.name,'OT'))
                filepath = [Brain_path,filesep,filename.name];
                files = dir(filepath);
                files = files(3:end);
                for kk = 1: length(files)
                    modality=files(kk);
                    if ~isempty(strfind(modality.name,'OT')) && isempty(strfind(modality.name,'mask')) && isempty(strfind(modality.name,'N4ITK'))
                        modality_path =  [Brain_path,filesep,filename.name,filesep,modality.name];
                        truth = mha_read_volume(modality_path);
                        info  = mha_read_header(modality_path);
                    end
                end 
            end  
             if ~isempty(strfind(filename.name,'T1')) && isempty(strfind(filename.name,'T1c'))
                filepath = [Brain_path,filesep,filename.name];
                files = dir(filepath);
                files = files(3:end);
                for kk = 1: length(files)
                    modality=files(kk);
                    if ~isempty(strfind(modality.name,'T1')) && isempty(strfind(modality.name,'mask')) && isempty(strfind(modality.name,'N4ITK'))
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
                    if ~isempty(strfind(modality.name,'T2')) && isempty(strfind(modality.name,'mask')) && isempty(strfind(modality.name,'N4ITK'))
                        modality_path =  [Brain_path,filesep,filename.name,filesep,modality.name];
                        T2 = mha_read_volume(modality_path);
                        info  = mha_read_header(modality_path);
                    end
                end 
            end  
            if ~isempty(strfind(filename.name,'Flair'))
                filepath = [Brain_path,filesep,filename.name];
                files = dir(filepath);
                files = files(3:end);
                for kk = 1: length(files)
                    modality=files(kk);
                    if ~isempty(strfind(modality.name,'Flair')) && isempty(strfind(modality.name,'mask')) && isempty(strfind(modality.name,'N4ITK'))
                        flair_name = modality.name;
                        modality_path =  [Brain_path,filesep,filename.name,filesep,modality.name];
                        FLAIR = mha_read_volume(modality_path);
                        info_flair  = mha_read_header(modality_path);
                    end
                end 
            end  
            if ~isempty(strfind(filename.name,'T1c'))
                filepath = [Brain_path,filesep,filename.name];
                files = dir(filepath);
                files = files(3:end);
                for kk = 1: length(files)
                    modality=files(kk);
                    if ~isempty(strfind(modality.name,'T1c')) && isempty(strfind(modality.name,'mask')) && isempty(strfind(modality.name,'N4ITK'))
                        modality_path =  [Brain_path,filesep,filename.name,filesep,modality.name];
                        T1C = mha_read_volume(modality_path);
                        info  = mha_read_header(modality_path);
                    end
                end 
            end  
            
  if compute_statistics ==0           
   truth= NaN;        
  end
end      
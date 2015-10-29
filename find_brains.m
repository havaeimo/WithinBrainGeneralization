brain_dir_ref = ['/home/local/USHERBROOKE/havm2701/data/Brats_train_2014/'];


brain_type_list_ref = dir(brain_dir_ref);
brain_type_list_ref = brain_type_list_ref(3:end);

 for subfoldercounter_ref=1:length(brain_type_list_ref)
     type_ref = brain_type_list_ref(subfoldercounter_ref).name;
     if type_ref == 'LG'
         continue;
     end
     brain_type_path_ref = [brain_dir_ref,type_ref];
     brain_list_ref = dir(brain_type_path_ref);
     brain_list_ref =brain_list_ref(3:end);
    
     for brain_counter_ref = 1:length(brain_list_ref)
         name = brain_list_ref(brain_counter_ref).name;
         brain_name_path_ref = [brain_type_path_ref,'/',name];
 
        if ~isdir(brain_name_path_ref)
            continue
        end        

%% load modalities
        [T1C_ref, T2_ref, FLAIR_ref, truth_ref, info_ref,flair_name_ref] = load_modalities_mac(brain_name_path_ref,0);
        truth_ref = zeros(size(T2_ref,1),size(T2_ref,2),size(T2_ref,3));
        [height_ref width_ref depth_ref] = size(T1C_ref);
        %T1C_ref = double(T1C_ref);
        %T2_ref = double(T2_ref);
        %FLAIR_ref = double(FLAIR_ref);
        % take backup of the modalities
        %%

                    brain_dir_qur = ['/home/local/USHERBROOKE/havm2701/data/mha_data_2013/Challenge_mha/'];


                    brain_type_list_qur = dir(brain_dir_qur);
                    brain_type_list_qur = brain_type_list_qur(3:end);

                     for subfoldercounter_qur=1:length(brain_type_list_qur)
                         type_qur = brain_type_list_qur(subfoldercounter_qur).name;
                          if type_qur == 'LG'
                             continue;
                         end
                         brain_type_path_qur = [brain_dir_qur,type_qur];
                         brain_list_qur = dir(brain_type_path_qur);
                         brain_list_qur =brain_list_qur(3:end);

                         for brain_counter_qur = 1:length(brain_list_qur)
                             name = brain_list_qur(brain_counter_qur).name;
                             brain_name_path_qur = [brain_type_path_qur,'/',name];

                            if ~isdir(brain_name_path_qur)
                                continue
                            end

                    %% load modalities
                            [T1C_qur, T2_qur, FLAIR_qur, truth_qur, info_qur,flair_name_qur] = load_modalities_mac(brain_name_path_qur,0);
                            truth_qur = zeros(size(T2_qur,1),size(T2_qur,2),size(T2_qur,3));
                            [height_qur width_qur depth_qur] = size(T1C_qur);
                            %T1C_qur = double(T1C_qur);
                            %T2_qur = double(T2_qur);
                            %FLAIR_qur = double(FLAIR_qur);
                            % take backup of the modalities



  %%
  if height_qur==height_ref && width_qur==width_ref && depth_qur==depth_ref
  
                check_T1C = T1C_qur - T1C_ref ; 
                if sum(check_T1C(:)) < 10

                    display(brain_name_path_ref)
                end
  end

                         end
                      end  
     end
 end
 
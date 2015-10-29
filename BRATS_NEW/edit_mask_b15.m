function MASK = edit_mask(T1C, T2, FLAIR,MASK)
[height,width,depth] = size(T1C);
display('****Edit the mask*****')
while 1
    %MASK_slice = imread('Mask_HG14_slice80.png');
    while 1
        fig2=figure();
         d=input('press "1" to select healthy region');
        if d==1
            sliceH = input('enter slice number for healthy (from image J):');
%             sliceH = depth - sliceH;
            modality = input('Enter modality to slelect from: ','s');
            view = input('Enter view to slelect from: ');
            imshow_modality(sliceH,modality,view);
            
            roi= roipoly;
            if view ~= 3 
                roi = flipud(roi)';
            elseif view ==3
                roi = roi';
            end
            MASK_sliceH = roi;
            MASK_sliceH = 10 * MASK_sliceH;
            imshow(roi);
            roi_index = find(roi>0);

            if view ==1
            MASK_sliceH =MASK(sliceH,:,:);
            MASK_sliceH(roi_index)= 10;
            MASK(sliceH,:,:) = MASK_sliceH;
            elseif view ==2
            MASK_sliceH =MASK(:,sliceH,:);
            MASK_sliceH(roi_index)= 10;
            MASK(:,sliceH,:) = MASK_sliceH;
            elseif view == 3
            MASK_sliceH =MASK(:,:,sliceH);
            MASK_sliceH(roi_index)= 10;
            MASK(:,:,sliceH) = MASK_sliceH;
            end
        else 
            break
        end                
    end
    while 1
        d=input('press "1" to select edema region');
        if d==1
            sliceE = input('enter slice number for edema (from image J):');
%             sliceE = depth - sliceE;
            modality = input('Enter modality to slelect from: ','s');
            view = input('Enter view to slelect from: ');
            imshow_modality(sliceE,modality,view);
             roi= roipoly;
            if view ~= 3 
                roi = flipud(roi)';
            elseif view ==3
                roi = roi';
            end
            imshow(roi);
            roi_index = find(roi>0);

            if view ==1
            MASK_sliceE =MASK(sliceE,:,:);
            MASK_sliceE(roi_index)= 2;
            MASK(sliceE,:,:) = MASK_sliceE;
            elseif view ==2
            MASK_sliceE =MASK(:,sliceE,:);
            MASK_sliceE(roi_index)= 2;
            MASK(:,sliceE,:) = MASK_sliceE;
            elseif view == 3
            MASK_sliceE =MASK(:,:,sliceE);
            MASK_sliceE(roi_index)= 2;
            MASK(:,:,sliceE) = MASK_sliceE;
            end
        else 
            break
        end                
    end


    dd=input('press "1" to continue selecting ROI');
     if dd~=1
          break
     end
end
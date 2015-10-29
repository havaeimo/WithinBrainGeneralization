function MASK = edit_mask_brats(MASK)
[height,width,depth] = size(MASK);
display('****Edit the mask*****')
while 1
    %MASK_slice = imread('Mask_HG14_slice80.png');
   
        fig2=figure();
         d=input('press label to select ');

            sliceH = input('enter slice number (from image J):');
%             sliceH = depth - sliceH;

            view = input('Enter view to slelect from: ');
            imshow_predict(MASK,sliceH,view);
            
            roi= roipoly;
            if view ~= 3 
                roi = flipud(roi)';
            elseif view ==3
                roi = roi';
            end
            MASK_sliceH = roi;
            MASK_sliceH = d * MASK_sliceH;
            imshow(roi);
            roi_index = find(roi>0);

            if view ==1
            MASK_sliceH =MASK(sliceH,:,:);
            MASK_sliceH(roi_index)= d;
            MASK(sliceH,:,:) = MASK_sliceH;
            elseif view ==2
            MASK_sliceH =MASK(:,sliceH,:);
            MASK_sliceH(roi_index)= d;
            MASK(:,sliceH,:) = MASK_sliceH;
            elseif view == 3
            MASK_sliceH =MASK(:,:,sliceH);
            MASK_sliceH(roi_index)= d;
            MASK(:,:,sliceH) = MASK_sliceH;
            end
              
    dd=input('press "1" to continue selecting ROI');
     if dd~=1
          break
     end
    
end
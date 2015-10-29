function imshow_predict(MASK, sliceE,view)

[height,width,depth] = size(MASK);
    
        if view ==3
            imshow(MASK(:,:,sliceE)',[])
        elseif view == 1
            imshow(flipud(reshape(MASK(sliceE,:,:),width,depth)'),[])
        elseif view == 2
             imshow(flipud(reshape(MASK(:,sliceE,:),height,depth)'),[])
        end

end


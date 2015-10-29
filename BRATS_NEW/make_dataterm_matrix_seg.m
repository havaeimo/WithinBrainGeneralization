function matrix = make_dataterm_matrix_seg(IDX,mask_idx,space,T1C,num_labels)

k = size(IDX,2);

labels = 0:num_labels-1;

matrix = zeros(length(labels),length(T1C(:)));
matrix(1,:) = k;

[h,w,d] = size(T1C);
for count = 1: length(IDX(:,1))

        ROW = uint32(space(2,count) * h);
        COL = uint32(space(3,count) * w);
        DEP = uint32(space(4,count) * d);
       % importance(ROW,COL,DEP) = min(D(count,:));
        index = (h*w)*(DEP-1)+(COL - 1)*h+ROW;

        M = mask_idx(IDX(count,:));
        M(M==10)=0;
        M(M==3)=1;
        M(M==4)=3;
   
        matrix(:,index) = hist(M,labels);
    %    end
end 




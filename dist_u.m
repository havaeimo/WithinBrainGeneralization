function d = dist_u(x,y)
  %d = zeros(length(ZJ),1);
  %for i = 1:length(ZJ)
    %y = ZJ(i,:);  
    sum = 0;
    sum = sum+ (x(1)-y(1)).^2;
    sum = sum+ (x(2)-y(2)).^2;
    sum = sum+ (x(3)-y(3)).^2;
    sum = sum+ (1 - gaussmf(abs(x(4)-y(4)),[0.5 0]));
    sum = sum+ (1 - gaussmf(abs(x(5)-y(5)),[0.5 0]));
    sum = sum+ (1 - gaussmf(abs(x(6)-y(6)),[0.5 0]));
    %d(i) = sqrt(sum);
    d = sqrt(sum);
    %end
    
function [label,spatial_info,prob] = read_libsvm_probability_line_2(line,classes)

a = strread(line, '%s', 'delimiter', ' ');
a = str2double(a);
label = a(1);
spatial_info = a(end-2:end);

prob = a(2:6);

prob(end) = [] ; % helathy class has label zero. so by removing it prob(ind) corrstponds to class ind







% 
% for i=1:length(a)-1
%     output(i) = str2double(a(i+1));
% end
% 
% output(end)=[];
% 
% [label x_str y_str z_str] = strread(line, '%s %s %s %s', 'delimiter', ' ');
% [ix x] = strtok(x_str,'0');
% [iy y] = strtok(y_str, '0');
% [iz z] = strtok(z_str, '0');
% x = str2double(x);
% y = str2double(y);
% z = str2double(z);
% label = str2double(label);
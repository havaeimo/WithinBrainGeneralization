function all_metrics = main()
% bmvc find parameters
%main_path = 'N:\User\Research\Results\BMVC_2014_find_parameters\Tree_lpboost\Newmasks_6Dfeature\boykov\Denoised';
main_path = 'N:\User\Research\Results\BMVC_2014_find_parameters\kNN_rbf\Newmasks_6Dfeature\boykov\Denoised\';

parameter_paths = dir(main_path);
parameters = parameter_paths(3:end);
nb_parameters = length(parameters);

% make a matrix size (number_parametercombinations) x (3) c [alpha,beta,avg_f1score,total_mean]
all_metrics = zeros(nb_parameters,4);

for parameter = 1:nb_parameters
    parameter_name = parameters(parameter).name;
    parameter_path = [main_path,parameter_name,'\'];
    
    [alpha , beta] = find_alpha_beta(parameter_name);
    [f1_score_average,total_mean_vector] = fscore_parameters(parameter_path);
    all_metrics(parameter,:) = [alpha, beta, f1_score_average,total_mean_vector];
end


end

function f1_avg = F_score(stat_measure)
    Pr = stat_measure(:,1);
        Re = stat_measure(:,2);
        sum = 0;beta = 0;
        for l = 1:4
            value= 2*(Pr(l).*Re(l))./(Pr(l)+Re(l));
           if isnan(value)
                beta = beta+1;
                continue
           end
           sum = sum+value;    
        end
        f1_avg=sum/4-beta;
           
end


function avg = total_avg(stat_measure)
    stat_measure = stat_measure(:);
    sum = 0;beta = 0;
    for l = 1:16
        value = stat_measure(l);
           if isnan(value)
                beta = beta+1;
                continue
           end
           sum = sum+value;    
    end
    avg = sum/16-beta;
end

function [alpha, beta] = find_alpha_beta(string)
    [c,a,b] = strread(string,'%s %f %f','delimiter','evaluationmatrixalpha_' );
    alpha = max(a);
    beta = max(b);
end

function [f1_score,total_mean] = fscore_parameters(path)
    brain_paths =  dir(path);
    brain_paths = brain_paths(3:end);
    nb_brains = length(brain_paths);
    
    f1_vector = zeros(1,nb_brains);
    total_mean_vector = zeros(1,nb_brains);
    for i= 1:nb_brains
        stat = load([path,brain_paths(i).name]);
        f1 = F_score(stat.evaluation_matrix);
        m = total_avg(stat.evaluation_matrix);
        f1_vector(i) = f1;
        total_mean_vector(i) = m;
    end
    f1_score = mean(f1_vector); % the f1_score average across all brains with parameters alpha beta
    total_mean = mean(total_mean_vector);
end


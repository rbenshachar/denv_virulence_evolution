function param_estimates = MLE_DENV_evolution_subset_data_innate()

load('low_params')

time_PI = params.time_PI_low;  
data_PI = params.data_PI_low; 
LOD_PI = params.LOD_PI_low; 
 
       %take out all Nans
    for i = 1:11
        index = isnan(data_PI(i,:));
        time_temp = time_PI(i,:);
        data_temp  = data_PI(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_PI(i).time = time_temp; 
        indiv_PI(i).data = data_temp; 
        indiv_PI(i).LOD =  LOD_PI(i,1);
        indiv_PI(i).c = zeros(1, length(indiv_PI(i).data)); 
             for j = 1:length(indiv_PI(i).data)
                        if data_temp(j) <= indiv_PI(i).LOD      
                             indiv_PI(i).c(j) = 0;  
                        else
                             indiv_PI(i).c(j)  = 1; 
                        end
            end
    end
 
time_SI = params.time_SI_low; 
data_SI = params.data_SI_low; 
LOD_SI = params.LOD_SI_low; 
 
  for i = 1:90
        index = isnan(data_SI(i,:));
        time_temp = time_SI(i,:);
        data_temp  = data_SI(i,:);
        data_temp(isnan(data_temp)) = [];
        time_temp(index) = [];
        indiv_SI(i).time = time_temp; 
        indiv_SI(i).data = data_temp; 
        indiv_SI(i).LOD =  LOD_SI(i,1);
        indiv_SI(i).c = zeros(1, length(indiv_SI(i).data)); 
             for j = 1:length(indiv_SI(i).data)
                        if data_temp(j) <= indiv_SI(i).LOD      
                             indiv_SI(i).c(j) = 0;  
                        else
                             indiv_SI(i).c(j)  = 1; 
                        end
            end
  end
  
parameters = [params.time_start, params.time_end, params.Xinit, params.Yinit, params.Vinit...
 params.Ninit, params.alpha, params.omega, params.d, 1e-10, 20,1e-4, params.deltaT,params.dT,params.Tinit,1e-3, 5.5, 0.15, 1];
  
load('param_estimates_full.mat');
a = (param_estimates(1) + param_estimates(4))/param_estimates(1);
b = param_estimates(5); 
vals = [a, b];

%qPI, qT, IP
init_params = cat(2, 0.1, .01, 4e-10);

[param_estimates,  f_val] = fminsearch(@(x)MLE_likelihood_subset_data_innate(x, parameters, indiv_PI, indiv_SI, vals), init_params); %, options);
save('param_estimates_low_innate', 'param_estimates')
save('fval_low_innate', 'f_val')


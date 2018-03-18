function param_estimates = MLE_DENV_evolution_full()

%estimate parameters for full model

load('params')
params.d = 0.07; 
params.dT = 0; 

time_PI = params.time_PI; 
data_PI = params.data_PI; 
LOD_PI = params.LOD_PI; 
 
    %take out all Nans
    for i = 1:30
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
  
 time_SI = params.time_SI; 
 data_SI = params.data_SI; 
 LOD_SI = params.LOD_SI;
 
 
  for i = 1:209
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
 
init_params = [5e-10, 4e-4, 4e-07, 1e-10,  0.15, 0.07];

[param_estimates, f_val] = fminsearch(@(x)MLE_likelihood_full(x, parameters, indiv_PI, indiv_SI), init_params); %, options);
save('param_estimates_full_test', 'param_estimates')
save('f_val_full_test', 'f_val')


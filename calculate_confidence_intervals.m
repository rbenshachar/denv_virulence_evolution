function CI = calculate_confidence_intervals(model)

% calculating the 95% confidence intervals (CI)
% we calculate CI as parameter values that result in change of +/-2 likelihood 

%model = 1 : model fit to full dataset
%model = 2 : model fit to subset of data
%model = 3 : innate IR model fit subset of data
%model = f : alternative T-cell model a = 10^6
%model = 5: alternative T-cell model a = 10^7

if model == 1 
    load('param_estimates_full.mat');
    load('f_val_full')
elseif model == 2
    load('low_params');
    load('param_estimates_low.mat');
    load('fval_low')
elseif model == 3
    load('low_params');
    load('param_estimates_low_innate.mat');
    load('fval_low_innate')
elseif model == 4
    load('low_params');
    load('param_estimates_low_alt_T_a_1e6.mat');
    load('fval_low_alt_T_a_1e6');
elseif model == 5
    load('low_params');
    load('param_estimates_low_alt_T_a_1e7.mat');
    load('fval_low_alt_T_a_1e7');
end

n_params = length(param_estimates);
RMS = f_val;

LB = zeros(length(n_params));
UB = zeros(length(n_params));

params_orig = param_estimates; 

if model == 1
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
  
else
    time_PI = cat(1, params.time_PI_1_low, params.time_PI_2_low, params.time_PI_3_low);  
    data_PI = cat(1, params.data_PI_1_low, params.data_PI_2_low, params.data_PI_3_low); 
    LOD_PI = cat(1, params.LOD_PI_1_low, params.LOD_PI_2_low, params.LOD_PI_3_low); 
 
    params.d = 0; 
 
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
 
    time_SI = cat(1, params.time_SI_1_low, params.time_SI_2_low, params.time_SI_3_low, params.time_SI_4_low); 
    data_SI = cat(1, params.data_SI_1_low, params.data_SI_2_low, params.data_SI_3_low, params.data_SI_4_low); 
    LOD_SI = cat(1, params.LOD_SI_1_low, params.LOD_SI_2_low, params.LOD_SI_3_low, params.LOD_SI_4_low); 
 
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
  
end

 parameters = [params.time_start, params.time_end, params.Xinit, params.Yinit, params.Vinit...
        params.Ninit, params.alpha, params.omega, params.d, 1e-10, 20,1e-4, params.deltaT,params.dT,params.Tinit,1e-3, 5.5, 0.15, 1];

for i = 1:n_params
    v = param_estimates; 
    [LB(i),  ~] = fminsearch(@(x)f1(i, v, x, parameters, indiv_PI, indiv_SI, RMS, model), v(i)); 
    v = param_estimates; 
    [UB(i),  ~] = fminsearch(@(x)f2(i, v, x, parameters, indiv_PI, indiv_SI, RMS, model), v(i)); 
end

CI = [LB, UB]; 
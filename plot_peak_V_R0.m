function [R0_omega, disease_omega, max_clear_omega] = plot_peak_V_R0(s, m, k)

%FIT TO SUBSET OF DATA

%s = serotype
%m = immune status
%k = number of simulations

load('params')
params.d = 0; 
params.dT = 0; 

load('param_estimates_full.mat');
params.cv = param_estimates(5); 

params.sigmae = 0.5; 

load('param_estimates_low.mat')
params.beta_PI = param_estimates(3);
params.q_PI = param_estimates(1); 
params.beta_SI = param_estimates(3)*1.27; 

params.duration_cut_off = log10(1500);  
        
R0_omega = zeros(k, 1);
disease_omega = R0_omega; 
max_clear_omega = R0_omega; 

if m == 1

    for i = 1:k 
        %omega
        params.time_end = 25; 
        params.beta = params.beta_PI; 
        params.q = params.q_PI;

        phenotype = params.omega; 
        params.beta = normrnd(params.beta_PI, params.beta_PI*params.cv);        
    
        [disease_omega(i,1), ~, R0_omega(i,1), ~,max_clear_omega(i,1)] = calculate_PI_fitness_disease(phenotype, s, params, 1/24);
  
        i
    end

elseif m == 2

    for i = 1:k 

        params.time_end = 50; 
        params.q = 0;  
        params.qT = param_estimates(2);  

        phenotype = params.omega; 
        params.beta = normrnd(params.beta_SI, params.beta_SI*params.cv);
 
    [R0_omega(i,1), disease_omega(i,1), max_clear_omega(i,1)] = calculate_SI_fitness_disease(phenotype,s, params, 1/24); 
    i
    end

end
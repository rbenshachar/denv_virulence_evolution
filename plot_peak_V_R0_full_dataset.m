function [R0_omega, disease_omega, max_clear_omega] = plot_peak_V_R0_full_dataset(s, m, k)

%s = serotype
%m = immune status
%k = number of simulations

load('params')

params.sigmae = 0.5; 
params.dT = 0; 
params.duration_cut_off = log10(1500);  

params.omega = 1e4; 
params.kappa = 5; 
params.Vinit = 10^(-3);
        
R0_omega = zeros(k, 1);
disease_omega = R0_omega; 
max_clear_omega = R0_omega; 

load('param_estimates_full.mat');
params.q_PI = param_estimates(2);
params.cv = param_estimates(5); 
params.d = param_estimates(6);

if m == 1
    for i = 1:k 
        %omega
        params.time_end = 25;
        params.q = params.q_PI;
        
        phenotype = params.omega; 
        beta = param_estimates(1); 
        params.beta = normrnd(beta, params.cv*beta);

        [disease_omega(i,1), ~, R0_omega(i,1), ~,max_clear_omega(i,1)] = calculate_PI_fitness_disease(phenotype,s, params, 1/24);   
    i
    end
else

    for i = 1:k
        %omega
        params.time_end = 40; 
        params.q = 0;  
        params.qT = param_estimates(3);  
        
        phenotype = params.omega; 
        beta =  param_estimates(1) + param_estimates(4); 
        params.beta = normrnd(beta, params.cv*beta);
        
        [R0_omega(i,1), disease_omega(i,1), max_clear_omega(i,1)]  = calculate_SI_fitness_disease(phenotype, s, params, 1/24);
        i
    end
end

end


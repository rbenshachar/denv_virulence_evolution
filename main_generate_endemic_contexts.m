function void = main_generate_endemic_contexts(void)

% this function generates outfiles of endemic contexts for either the no
% cross-immunity case, the classical cross-immunity case, or the clinical
% cross-immunity case. The no cross-immunity case is just implemented by setting delta (the waning of immunity from T to S class) 
% to be a very high value, corresponding to a very short duration of cross-protection

% 1 = clinical cross-protection
% 0 = classical cross-protection

% delta = rate of waning heterologous immunity

clear all; close all; clc; tic

infile_TP = 'figure4_infile'; load(infile_TP);
  
% start with some infile to provide initial conditions:
dengue_infile = 'dengue_init_file_DENV1'; 

% Plot infecteds for this infile:
%PlotInfecteds(dengue_infile); % the dynamics in this infile are at equilibrium

% now, generate endemic contexts for a range of f and omega values

t_sim_yrs = 20; 
for omega_cell = 65:170

    for f = 1:0.5:3
    
        R0_PI = f*TP_primary(1,omega_cell);
        R0_SI = f*TP_secondary(1,omega_cell);

        % simulates for t_sim_yrs years to generate new endemic context at equilibrium
        dengue_outfile = strcat('dengue_init_file_DENV1_f', int2str(f*100), '_w', int2str(omega_vals(omega_cell)), '_noCrossProtection')
        RunDengueAgeSimulations(R0_PI, R0_SI, -t_sim_yrs, 0, dengue_outfile, dengue_infile, 0, 1/365); % one day cross-protection (NO CROSS-PROTECTION)
                
        dengue_outfile = strcat('dengue_init_file_DENV1_f', int2str(f*100), '_w', int2str(omega_vals(omega_cell)), '_classicalProtection')
        RunDengueAgeSimulations(R0_PI, R0_SI, -t_sim_yrs, 0, dengue_outfile, dengue_infile, 0, 2); % two years cross-protection (CLASSICAL CROSS-PROTECTION)
        
        dengue_outfile = strcat('dengue_init_file_DENV1_f', int2str(f*100), '_w', int2str(omega_vals(omega_cell)), '_clinicalProtection')
        RunDengueAgeSimulations(R0_PI, R0_SI, -t_sim_yrs, 0, dengue_outfile, dengue_infile, 1, 2); % two years cross-protection (CLINICAL CROSS-PROTECTION)
        
    end
    
end

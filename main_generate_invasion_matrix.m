function void = main_generate_invasion_matrix(void)

% this code generates a PIP matrix for a given transmission factor f, going
% from omega value in omega_cell_start cell number to omega_cell_end cell
% number

clear all; close all; clc; tic

infile_TP = 'figure4_infile'; load(infile_TP); 

f = 2; omega_cell_start = 65; omega_cell_end = 170;

omega_cell_list = omega_cell_start:omega_cell_end;
PIP_matrix = NaN*ones(length(omega_cell_list));

outfile = strcat('PIP_noCrossProtection_f', int2str(f*100));
%outfile = strcat('PIP_classicalProtection_f', int2str(f*100));
%outfile = strcat('PIP_clinicalProtection_f', int2str(f*100));

omega_cntr_res = 1;
for omega_cell_res = omega_cell_list
    
    dengue_infile = strcat('dengue_init_file_DENV1_f', int2str(f*100), '_w', int2str(omega_vals(omega_cell_res)), '_noCrossProtection')
    %dengue_infile = strcat('dengue_init_file_DENV1_f', int2str(f*100), '_w', int2str(omega_vals(omega_cell_res)), '_classicalProtection')
    %dengue_infile = strcat('dengue_init_file_DENV1_f', int2str(f*100), '_w', int2str(omega_vals(omega_cell_res)), '_clinicalProtection')
    
    R0_PI_res = f*TP_primary(1,omega_cell_res);
    R0_SI_res = f*TP_secondary(1,omega_cell_res);
            
    cumulativeI_res = RunDengueAgeSimulations_invasion(R0_PI_res, R0_SI_res, 0, 1, dengue_infile);
    
    omega_cntr_inv = 1;
    for omega_cell_invade = omega_cell_list
        
        if omega_cntr_res == omega_cntr_inv
            PIP_matrix(omega_cntr_inv,omega_cntr_res) = 1;
        else
            R0_PI_in = f*TP_primary(1,omega_cell_invade);
            R0_SI_in = f*TP_secondary(1,omega_cell_invade);
            
            cumulativeI_inv = RunDengueAgeSimulations_invasion(R0_PI_in, R0_SI_in, 0, 1, dengue_infile);
           
            % because of slight fluctuations in the number of susceptibles
            % over a year (due to birth/death accounting), it is simpler to
            % fix S and T values at the endemic equibrium at the very
            % beginning of a year, and then compare cumulative number of
            % infecteds of residents against the invading phenotype. A
            % phenotype that would be able to invade when rare would result
            % in a higher number of cumulative infections relative to the
            % resident strain.
            if cumulativeI_inv > cumulativeI_res % then can_increase_when_rare
                PIP_matrix(omega_cntr_inv,omega_cntr_res) = 1;
            else
                PIP_matrix(omega_cntr_inv,omega_cntr_res) = 0;
            end
                
        end
        omega_cntr_inv = omega_cntr_inv+1;
    end
    
    omega_cntr_res = omega_cntr_res + 1;
    PIP_matrix
    
    save(outfile, 'omega_vals', 'omega_cell_list', 'PIP_matrix');
end
save(outfile, 'omega_vals', 'omega_cell_list', 'PIP_matrix');
figure; pcolor(omega_vals(omega_cell_list), omega_vals(omega_cell_list), PIP_matrix); 
xlabel('\omega of resident strain'); ylabel('\omega of invading strain'); title(strcat('f', int2str(f*100)));

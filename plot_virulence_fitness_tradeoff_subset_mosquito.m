function plot_virulence_fitness_tradeoff_subset_mosquito()

%%%%%%%%%%%%%%%%%
%FIGURE S8 
%plotting effect of omega on mosquito

close all;

load('params');

colors = colormap(cbrewer('qual', 'Dark2', 5));

phenotype = 1e4*(.01:.01:4); 

eps0 = -10;
eps1 = 1e-3; 
phm =  1./(1 + exp(-eps1.*phenotype - eps0));

params.sigmae = 0; 
params.d = 0;  
params.duration_cut_off = log10(1500); 

params.omega = 1e4; 
params.kappa = 5; 
params.Vinit = 10^(-3);

load('param_estimates_full.mat')
f = (param_estimates(1) + param_estimates(4))/param_estimates(1); 

load('param_estimates_low.mat')
params.time_end = 300; 
params.beta = param_estimates(3); 
params.q = param_estimates(1); 

[~, ~, R0_PI_mean_1, ~, ~] = calculate_PI_fitness_disease(phenotype, 1, params, .01);

% % % % % % % % % SI
params.time_end = 500; %210; 
params.beta =  (param_estimates(3))*f;
params.q = 0;  
params.qT = param_estimates(2); 
    
[R0_SI_mean_1, ~, ~] = calculate_SI_fitness_disease(phenotype, 1, params, .01);

x = 1; y = 2; 

subplot(x,y,1)
plot(phenotype, phm, 'LineWidth', 2, 'Color', 'k')
xlabel('viral production rate (genome copies/cell/day)'); ylabel('p_{hm}')
text(0.1,0.9,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
set(gca, 'FontSize', 14)

subplot(x,y,2)
plot(phenotype, R0_PI_mean_1.*phm','Color', colors(1,:), 'LineWidth', 2)
hold on; 
plot(phenotype, R0_SI_mean_1.*phm','Color', colors(2,:), 'LineWidth', 2)
xlabel('viral production rate (genome copies/cell/day)'); ylabel('component of R_0')
text(0.1,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
set(gca, 'FontSize', 14)

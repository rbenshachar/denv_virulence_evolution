function void = plot_Nguyen_fit(void)

%FIGURE S7
%Supplemental figure - fits to subset with different Nguyen fits
close all; 

colors = colormap(cbrewer('qual', 'Dark2', 5));

[final_params_1, ~] = fit_DENV_viral_load(1);

V = -2:.1:12; 

mu = final_params_1(1); 
sigma = final_params_1(2); 

q_1 =  1./(1 + exp(-sigma.*V - mu));
 
[final_params_2, ~] = fit_DENV_viral_load(2);
mu = final_params_2(1); 
sigma =  final_params_2(2); 

q_2 =  1./(1 + exp(-sigma.*V - mu));

[final_params_3, ~] = fit_DENV_viral_load(3);

mu = final_params_3(1); 
sigma = final_params_3(2); 

q_3 =  1./(1 + exp(-sigma.*V - mu));

[final_params_4, ~] = fit_DENV_viral_load(4);

mu = final_params_4(1);
sigma = final_params_4(2);

q_4 =  1./(1 + exp(-sigma.*V - mu));

x = 1; y = 3; 

subplot(x,y,1)
plot(V, q_1, 'Color', colors(1,:), 'LineWidth', 2); 
hold on; 
plot(V, q_2, 'Color', colors(2,:), 'LineWidth', 2); 
hold on; 
plot(V, q_3, 'Color', colors(3,:), 'LineWidth', 2); 
hold on; 
plot(V, q_4, 'Color', colors(4,:), 'LineWidth', 2); 
text(0.1,0.9,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
xlim([-2, 12])
xlabel('viral load (log genome copies per ml)')
ylabel('Proportion of blood-fed mosquitoes with DENV infection')
set(gca, 'FontSize', 14);

load('params');
colors = colormap(cbrewer('qual', 'Dark2', 5));

params.sigmae = 0; 
params.d = 0; 
params.duration_cut_off = log10(1500); 

load('param_estimates_full.mat')
f = (param_estimates(1) + param_estimates(4))/param_estimates(1); 

load('param_estimates_low.mat')

phenotype = 1e4*(.01:.005:6); 

%PI
params.time_end = 300; 
params.beta = param_estimates(3); 
params.q = param_estimates(1); 
[disease_PI_mean_1, ~, TP_PI_mean_1, ~, ~] = calculate_PI_fitness_disease(phenotype, 1, params, .1);
params.time_end = 1000; 
[disease_PI_mean_2, ~, TP_PI_mean_2, ~, ~] = calculate_PI_fitness_disease(phenotype, 2, params, .1);
params.time_end = 300; 
[disease_PI_mean_3, ~, TP_PI_mean_3, ~, ~] = calculate_PI_fitness_disease(phenotype, 3, params, .1);
params.time_end = 1000; 
[disease_PI_mean_4, ~, TP_PI_mean_4, ~, ~] = calculate_PI_fitness_disease(phenotype, 4, params, .1);

% % % % % % % % % SI
params.time_end = 500; 
params.beta =  (param_estimates(3))*f;
params.q = 0;  
params.qT = param_estimates(2); 
[TP_SI_mean_1, disease_SI_mean_1, ~] = calculate_SI_fitness_disease(phenotype, 1, params, .1);
params.time_end = 1000; 
[TP_SI_mean_2, disease_SI_mean_2, ~] = calculate_SI_fitness_disease(phenotype, 2, params, .1);
params.time_end = 500; 
[TP_SI_mean_3, disease_SI_mean_3, ~] = calculate_SI_fitness_disease(phenotype, 3, params, .1);
params.time_end = 1500; 
[TP_SI_mean_4, disease_SI_mean_4, ~] = calculate_SI_fitness_disease(phenotype, 4, params, .1);

subplot(x,y,2)
plot(phenotype, TP_PI_mean_1, 'Color', colors(1,:), 'LineWidth', 2)
hold on; 
plot(phenotype, TP_SI_mean_1, 'Color', colors(1,:), 'LineStyle', '--', 'LineWidth', 2);
hold on; 
plot(phenotype, TP_PI_mean_2, 'Color', colors(2,:), 'LineWidth', 2)
hold on; 
plot(phenotype, TP_SI_mean_2, 'Color', colors(2,:), 'LineStyle', '--', 'LineWidth', 2);
hold on; 
plot(phenotype, TP_PI_mean_3, 'Color', colors(3,:), 'LineWidth', 2)
hold on; 
plot(phenotype, TP_SI_mean_3, 'Color', colors(3,:), 'LineStyle', '--', 'LineWidth', 2);
hold on; 
plot(phenotype, TP_PI_mean_4, 'Color', colors(4,:), 'LineWidth', 2)
hold on; 
plot(phenotype, TP_SI_mean_4, 'Color', colors(4,:), 'LineStyle', '--', 'LineWidth', 2);
text(0.1,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
xlabel('viral production rate (genome copies per cell per day)'); ylabel('transmission potential')
set(gca, 'FontSize', 14)

subplot(x,y,3)
plot(disease_PI_mean_1, TP_PI_mean_1, 'Color', colors(1,:), 'LineWidth', 2)
hold on; 
plot(disease_SI_mean_1, TP_SI_mean_1, 'Color', colors(1,:), 'LineStyle', '--', 'LineWidth', 2);
hold on; 
plot(disease_PI_mean_2, TP_PI_mean_2, 'Color', colors(2,:), 'LineWidth', 2)
hold on; 
plot(disease_SI_mean_2, TP_SI_mean_2, 'Color', colors(2,:), 'LineStyle', '--', 'LineWidth', 2);
hold on; 
plot(disease_PI_mean_3, TP_PI_mean_3, 'Color', colors(3,:), 'LineWidth', 2)
hold on; 
plot(disease_SI_mean_3, TP_SI_mean_3, 'Color', colors(3,:), 'LineStyle', '--', 'LineWidth', 2);
hold on; 
plot(disease_PI_mean_4, TP_PI_mean_4, 'Color', colors(4,:), 'LineWidth', 2)
hold on; 
plot(disease_SI_mean_4, TP_SI_mean_4, 'Color', colors(4,:), 'LineStyle', '--', 'LineWidth', 2);
text(0.1,0.9,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
xlabel('peak viral load (log genome copies per ml)'); ylabel('transmission potential')
xlim([4, 10])
set(gca, 'FontSize', 14)

end


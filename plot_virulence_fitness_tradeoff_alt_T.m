function void = plot_virulence_fitness_tradeoff_alt_T()
% FIGURE S5
% SATURATING T-CELL FORMULATION

close all;

load('params');

colors = colormap(cbrewer('qual', 'Dark2', 5));

params.sigmae = 0; 
params.d = 0; 
params.dT = 0; 
params.duration_cut_off = log10(1500); 

params.omega = 1e4; 
params.kappa = 5; 
params.Vinit = 10^(-3);
dt = 0.01; 

phenotype = 1e4.*(.01:.005:1); 

load('param_estimates_full.mat')
f = (param_estimates(1) + param_estimates(4))/param_estimates(1); 

load('param_estimates_low.mat')

params.time_end = 210;  
params.beta =  (param_estimates(3))*f;
params.q = 0;  
params.qT = param_estimates(2); 

[T, Y] = ode45(@(t,y)SI(t, y, params),params.time_start:dt:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit, params.Tinit]);   
V = log10(Y(:,3)); 
params.time_end = 210;
[TP_SI_mean, disease_SI_mean, ~] = calculate_SI_fitness_disease(phenotype, 1, params, .01);

params.a = 1e6; 
load('param_estimates_low_alt_T_a_1e6')
params.qT = param_estimates(3); 
params.beta = f*param_estimates(1); 
[T, Y3] = ode45(@(t,y)SI_alt(t, y, params),params.time_start:dt:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit, params.Tinit]);   
V3 = log10(Y3(:,3)); 
params.time_end = 210; 
[TP_SI_mean_1, disease_SI_mean_1, ~] = calculate_SI_fitness_disease_alt(phenotype, 1, params, .01);

params.a = 1e7; 
load('param_estimates_low_alt_T_a_1e7')
params.qT = param_estimates(3); 
params.beta = f*param_estimates(1);
[T, Y4] = ode45(@(t,y)SI_alt(t, y, params),params.time_start:dt:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit, params.Tinit]);   
V4 = log10(Y4(:,3)); 
params.time_end = 210; 
[TP_SI_mean_2, disease_SI_mean_2, ~] = calculate_SI_fitness_disease_alt(phenotype, 1, params, .01);

x = 1; y = 3; 
subplot(x,y,1)
plot(T, V, 'LineWidth', 2); 
hold on; 
plot(T, V3, 'LineWidth', 2, 'LineStyle', ':');
hold on; 
plot(T, V4, 'LineWidth', 2, 'LineStyle', '--');
xlabel('time since infection (days)');
ylabel('viral load (log genome copies per ml)')
text(0.1,0.9,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
ylim([-4, 10])
set(gca, 'FontSize', 14); 

subplot(x,y,2)
plot(phenotype, disease_SI_mean, 'LineWidth',2)
hold on
plot(phenotype, disease_SI_mean_1, 'LineWidth',2, 'LineStyle', ':')
hold on; 
plot(phenotype, disease_SI_mean_2, 'LineWidth',2, 'LineStyle', '--')
xlabel('viral production rate (copies/cell/day)'); ylabel('peak viral load (log genome copies per ml)')
text(0.1,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
set(gca, 'FontSize', 14); 
    
subplot(x,y,3)
plot(phenotype, TP_SI_mean,'LineWidth',2)
hold on
plot(phenotype, TP_SI_mean_1, 'LineWidth',2, 'LineStyle', ':')
hold on; 
plot(phenotype, TP_SI_mean_2, 'LineWidth',2, 'LineStyle', '--')
xlabel('viral production rate (genome copies per cell per day)'); ylabel('transmission potential')
text(0.1,0.9,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
set(gca, 'FontSize', 14); 
    

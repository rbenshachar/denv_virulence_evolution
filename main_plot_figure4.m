function void = main_plot_figure4(void)

%code to plot fig. 4 in manuscript
clear all; close all; clc;

infile = 'PIP_1serotype_r_TP_DENV1_FINAL_smooth'; load(infile); 

subplot(1,3,1); 
pcolor(omega_vals, omega_vals, PIP_matrix);
xlabel('\omega of resident strain'); ylabel('\omega of invading strain');
set(gca, 'FontSize', 14); 
opt_omega_oneSerotype = omega_vals(211)
text(0.9,0.95,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)

f_vals = 0.5:0.5:3;
subplot(1,3,2); 
plot(f_vals, opt_omega_oneSerotype*ones(size(f_vals)), 'b', 'LineWidth',2);
hold on;
plot(f_vals, opt_omega_oneSerotype*ones(size(f_vals)), 'b.', 'LineWidth',2); 
plot(2, opt_omega_oneSerotype, 'bo'); 
axis([0.25 3.25 min(omega_vals) max(omega_vals)]);
xlabel('transmission intensity f');
ylabel('optimal \omega');

opt_omega_twoSerotypes_noCrossImmunity_cells = [159 98 78 74 73 73]; 
opt_omega_twoSerotypes_noCrossImmunity = omega_vals(opt_omega_twoSerotypes_noCrossImmunity_cells);
plot(f_vals, opt_omega_twoSerotypes_noCrossImmunity, 'r', 'LineWidth',2); hold on;
plot(f_vals, opt_omega_twoSerotypes_noCrossImmunity, 'r.', 'LineWidth',2);
plot(2, opt_omega_twoSerotypes_noCrossImmunity(4), 'ro', 'LineWidth',2); 
opt_omega_twoSerotypes_noCross_f2 = opt_omega_twoSerotypes_noCrossImmunity(4)

opt_omega_twoSerotypes_classicalCrossImmunity_cells = [162 101 79 74 72 71]; 
opt_omega_twoSerotypes_classicalCrossImmunity = omega_vals(opt_omega_twoSerotypes_classicalCrossImmunity_cells);
plot(f_vals, opt_omega_twoSerotypes_classicalCrossImmunity, 'g', 'LineWidth',2); hold on;
plot(f_vals, opt_omega_twoSerotypes_classicalCrossImmunity, 'g.', 'LineWidth',2);

opt_omega_twoSerotypes_clinicalCrossImmunity_cells = [165 116 104 105 109 114];
opt_omega_twoSerotypes_clinicalCrossImmunity = omega_vals(opt_omega_twoSerotypes_clinicalCrossImmunity_cells);
plot(f_vals, opt_omega_twoSerotypes_clinicalCrossImmunity, 'm', 'LineWidth',2); hold on;
plot(f_vals, opt_omega_twoSerotypes_clinicalCrossImmunity, 'm.', 'LineWidth',2);
text(0.9,0.95,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
set(gca, 'FontSize', 14); 

load('params');

colors = colormap(cbrewer('qual', 'Dark2', 5));

params.sigmae = 0; 
params.d = 0; 
params.dT = 0; %0.1; 
params.duration_cut_off = log10(1500); 
params.delta = 0.3; 

load('param_estimates_full.mat');
params.cv = param_estimates(5);
f = (param_estimates(1) + param_estimates(4))/param_estimates(1); 

load('param_estimates_low.mat');

params.beta_PI = param_estimates(3);
params.beta_SI = f*(param_estimates(3)); 
dt = 0.1; 

params.q_PI = param_estimates(1); 
params.qT = (param_estimates(2)); 

%one circulating serotype
params.beta = params.beta_PI;
params.q = params.q_PI; 
params.omega = 22600; 
params.time_end = 50; 
[~, Y] = ode45(@(t,y)PI(t, y, params),params.time_start:dt:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);   
V = log10(Y(:,3));  

params.q = 0; 
params.beta = params.beta_SI;
[~, Y2] = ode45(@(t,y)SI(t, y, params),params.time_start:dt:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);   
V2 = log10(Y2(:,3));  

%two circulating serotypes
params.beta = params.beta_PI;
params.q = params.q_PI; 
params.omega = 14380; 
params.time_end = 50; 
[~, Y3] = ode45(@(t,y)PI(t, y, params),params.time_start:dt:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);   
V3 = log10(Y3(:,3));  

params.q = 0; 
params.beta = params.beta_SI;
[T, Y4] = ode45(@(t,y)SI(t, y, params),params.time_start:dt:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);   
V4 = log10(Y4(:,3));   

subplot(1,3,3)
plot(T, V, 'Color', colors(1,:),'LineWidth',2)
hold on; 
plot(T, V2, 'Color', colors(2,:),'LineWidth',2)
hold on; 
plot(T, V3, 'Color', colors(1,:),'LineStyle', '--','LineWidth',2)
hold on; 
plot(T, V4, 'Color', colors(2,:),'LineStyle', '--','LineWidth',2);
ylim([-4, 11]); 
xlabel('time since infection (days)')
ylabel('viral load (log genome copies/ml)')
text(0.9,0.95,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
set(gca, 'FontSize', 14); ylim([-4, 10])


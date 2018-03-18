function void = plot_viral_load_dynamics() 

%FIGURE S2
%SUPPLEMENTAL FIGURE OF TRADE-OFF DYNAMICS

close all; 
load('params');

params.sigmae = 0; 
params.d = 0; 
params.duration_cut_off = log10(1500); 
params.delta = 0.3; 

load('param_estimates_low.mat');

params.beta_PI = param_estimates(3);
dt = 0.1; 

params.q_PI = param_estimates(1); 
params.qT = (param_estimates(2)); 

phenotype = 1e4.*(1.5:.01:2.5);
params.time_end = 50; 
params.beta = params.beta_PI;
params.q = params.q_PI; 
[~, ~, R0_PI, ~, ~] = calculate_PI_fitness_disease(phenotype, 1, params, 0.1);

omega_PI = phenotype(R0_PI == max(R0_PI));

%primary
params.beta = params.beta_PI;
params.q = params.q_PI; 
params.omega = omega_PI(1); 
omega_PI(1)
params.time_end = 15; 
[~, Y] = ode45(@(t,y)PI_with_decay(t, y, params),params.time_start:dt:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);   
V = log10(Y(:,3)); 

params.omega = omega_PI(1) - 1e4.*1;
[~, Y] = ode45(@(t,y)PI_with_decay(t, y, params),params.time_start:dt:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);   
V_l = log10(Y(:,3));  

params.omega = omega_PI(1) + 1e4.*1; 
[T, Y] = ode45(@(t,y)PI_with_decay(t, y, params),params.time_start:dt:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);   
V_h = log10(Y(:,3));  

params.time_end = 15; 
TP = calculate_TP_over_time(omega_PI(1), 2, params, 1);
TP_l = calculate_TP_over_time(omega_PI(1) - 1e4.*1, 1, params, 1);
TP_h = calculate_TP_over_time(omega_PI(1) + 1e4.*1, 1, params, 1);

x = 1; y = 2; 

subplot(x,y,1)
plot(T, V_l, 'LineWidth',2)
hold on; 
plot(T, V, 'LineWidth', 2);
hold on; 
plot(T, V_h, 'LineWidth', 2); 
ylim([-4, 11]); 
xlabel('time since infection (days)')
ylabel('viral load (log genome copies per ml)')
text(0.9,0.95,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
set(gca, 'FontSize', 14); ylim([-4, 10])

t = params.time_start:.05:15; 
subplot(x,y,2)
plot(t,TP_l, 'LineWidth',2)
hold on; 
plot(t,TP, 'LineWidth',2)
hold on; 
plot(t,TP_h, 'LineWidth',2)
xlabel('time since infection (days)')
ylabel('transmission potential')
text(0.9,0.95,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
set(gca, 'FontSize', 14); ylim([0,5])



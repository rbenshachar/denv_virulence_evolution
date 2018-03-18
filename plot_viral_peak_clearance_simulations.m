function [p1, r1, p2, r2, p3, r3, p4, r4] = plot_viral_peak_clearance_simulations(k) 

%FIGURE 2

close all; 
colors = colormap(cbrewer('qual', 'Paired', 10));

x = 2; y = 2; 
    
load('params')

params.d = 0.07; 

load('param_estimates_full.mat');

params.beta_PI = param_estimates(1); 
params.q_PI = param_estimates(2); 
params.cv = param_estimates(5);  

params.beta_SI = param_estimates(1) + param_estimates(4);
params.qT = param_estimates(3); 
params.q_SI = 0; 

params.duration_cut_off = log10(1500);  
params.omega = 1e4; 
params.dT = 0; 
params.Vinit = 1e-3; 
params.kappa = 5; 

parameters = [params.time_start, params.time_end, params.Xinit, params.Yinit, params.Vinit...
params.Ninit, params.alpha, params.omega, params.d, 1e-10, 20, 1e-3, params.deltaT,params.dT,params.Tinit,1e-3, 5.5, 0, 1];


%primary infection
params.time_end = 15; 
params.beta = params.beta_PI; 
params.q = params.q_PI;

[T_PI, Y_PI] = ode45(@(t,y)PI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);  

params.beta = params.beta_PI - 2*params.cv*params.beta_PI; 
[~, Y2_PI] = ode45(@(t,y)PI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);  

params.beta = params.beta_PI + 2*params.cv*params.beta_PI; 
[~, Y3_PI] = ode45(@(t,y)PI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);  

params.time_end = 20; 
params.beta = params.beta_SI;  
params.q = params.q_SI;  

[T_SI, Y_SI] = ode45(@(t,y)SI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);  
V = log10(Y_SI(:,3));

params.beta = params.beta_SI - 2*params.cv*params.beta_SI; 
[~, Y2_SI] = ode45(@(t,y)SI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);  

params.beta = params.beta_SI + 2*params.cv*params.beta_SI; 
[~, Y3_SI] = ode45(@(t,y)SI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);  

%finding IP_js
params.IP = 5.9;
params.sigma = 0.2;

p = [0.001, 0.999];
crit = logninv(p, log(params.IP), params.sigma);
LB = round(crit(1), 2); 
UB = round(crit(2), 2); 
IP = LB:.01:UB; 

%PI
parameters(10) = params.beta_PI; %beta
parameters(11) = 5; %kappa
parameters(12) = params.q_PI; 
parameters(13) = 0; %deltaT
parameters(16) = 0; %qT
parameters(5) = params.Vinit; 

data = params.data_PI; 
time = params.time_PI; 
temp_LOD = params.LOD_PI; 

IP_val_PI = finding_IPs(parameters, data, time, temp_LOD, params, IP);

subplot(x,y,1)
for j = 1:size(data, 1) 
   plot(params.time_PI(j,:) + IP_val_PI(j), params.data_PI(j,:),'k-.', 'LineWidth',1);
   hold on; 
end

plot(0:15, log10(1500).*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','--'); 
hold on;
plot(0:15, log10(15000).*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','--');
hold on; 
plot(T_PI, log10(Y_PI(:,3)), 'Color', colors(8,:), 'LineWidth', 2);
hold on;
plot(T_PI, log10(Y2_PI(:,3)), 'Color', colors(8,:), 'LineWidth', 2, 'LineStyle', '--');
hold on; 
plot(T_PI, log10(Y3_PI(:,3)), 'Color', colors(8,:), 'LineWidth', 2, 'LineStyle', ':');
text(0.8,0.95,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 16)
ylabel('viral load (log genome copies per ml)'); title('1^° infection'); ylim([-2, 12]);
xlabel('time since infection (days)')
set(gca, 'FontSize', 14)
xlim([0, 15]);

%SI
%finding IP_js
params.IP = 5.9;
params.sigma = 0.2;

p = [0.001, 0.999];
crit = logninv(p, log(params.IP), params.sigma);
LB = round(crit(1), 2); 
UB = round(crit(2), 2); 
IP = LB:.01:UB; 

parameters(10) = params.beta_SI; %beta
parameters(11) = 5; %kappa
parameters(12) = params.q_SI; %q

parameters(13) = 1e-6; %deltaT
parameters(16) = params.qT; %qT

data = params.data_SI; 
time = params.time_SI; 
temp_LOD = params.LOD_SI; 

IP_val_SI = finding_IPs(parameters, data, time, temp_LOD, params, IP);

subplot(x,y,2)
for j = 1:size(data, 1)
   plot(params.time_SI(j,:) + IP_val_SI(j), params.data_SI(j,:),'k-.', 'LineWidth',1);
   hold on; 
end
plot(0:15, log10(1500).*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','--'); 
hold on;
plot(0:15, log10(15000).*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','--'); 
hold on;  
plot(T_SI, log10(Y_SI(:,3)), 'Color', colors(10,:), 'LineWidth', 2);
hold on;
plot(T_SI, log10(Y2_SI(:,3)), 'Color', colors(10,:), 'LineWidth', 2, 'LineStyle','--' );
hold on; 
plot(T_SI, log10(Y3_SI(:,3)), 'Color', colors(10,:), 'LineWidth', 2, 'LineStyle',':' );
text(0.8,0.95,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 16)
title('2^° infection');  ylim([-2, 12]); xlim([0, 15]);
xlabel('time since infection (days)') 
set(gca, 'FontSize', 14); 
 
% %simulated relationships
% %full dataset
%Primary infections
   [~, peakV, viral_clearance] = plot_peak_V_R0_full_dataset(1,1, k);
   f1 = fitlm(peakV',viral_clearance', 'linear');
   a1 = f1.Coefficients.Estimate; 
   p1 = f1.Coefficients.pValue(2);
   r1 = f1.Coefficients.Estimate(2);    

   subplot(x,y,3)
   plot(peakV, viral_clearance , '.', 'Color', colors(1,:), 'MarkerSize', 15, 'MarkerFaceColor','none'); 
   hold on; 
 

%secondary infections
[~,peakV, viral_clearance] = plot_peak_V_R0_full_dataset(1, 2, k);
f2 = fitlm(peakV',viral_clearance', 'linear');
a2 = f2.Coefficients.Estimate; 
p2 = f2.Coefficients.pValue(2);
r2 = f2.Coefficients.Estimate(2);    
   
v = 3:14;
subplot(x,y,4)
plot(peakV, viral_clearance , '.', 'Color', colors(1,:), 'MarkerSize', 15, 'MarkerFaceColor','none'); 
hold on; 

%subset of the data
%Primary infections
[~, peakV, viral_clearance] = plot_peak_V_R0(1,1, k);        
subplot(x,y,3)
plot(peakV, viral_clearance , 'o', 'Color', colors(2,:), 'MarkerSize', 5, 'MarkerFaceColor','none'); 
f3 = fitlm(peakV',viral_clearance', 'linear');
a3 = f3.Coefficients.Estimate; 
p3 = f3.Coefficients.pValue(2);
r3 = f3.Coefficients.Estimate(2);  
hold on; 
plot(v, a3(2).*v + a3(1), 'k-.','LineWidth', 2);  

v = 3:14;
plot(v, a1(2).*v + a1(1), 'k','LineWidth', 2); 
text(0.8,0.9,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 16)
xlabel('peak viral load (log genome copies per ml)'); ylabel('maximum daily viral clearance rate (log genome copies per ml per day)') 
set(gca, 'FontSize', 14); xlim([5, 12]); ylim([1, 8])

%secondary infections
[~, peakV, viral_clearance] = plot_peak_V_R0(1,2, k);        
subplot(x,y,4)
f4 = fitlm(peakV',viral_clearance', 'linear');
a4 = f4.Coefficients.Estimate; 
p4 = f4.Coefficients.pValue(2);
r4 = f4.Coefficients.Estimate(2);  
plot(peakV, viral_clearance , 'o', 'Color', colors(2,:), 'MarkerSize', 5, 'MarkerFaceColor','none'); 
hold on; 
plot(v, a4(2).*v + a4(1), 'k-.','LineWidth', 2);
hold on;
plot(v, a2(2).*v + a2(1), 'k','LineWidth', 2); 
text(0.8,0.9,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 16)
xlabel('peak viral load (log genome copies per ml per day)'); 
set(gca, 'FontSize', 14); xlim([5, 12]); ylim([1, 8])

end

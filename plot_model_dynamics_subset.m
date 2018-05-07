function void = plot_model_dynamics_subset(void)

%FIGURE S1

close all;

colors = colormap(cbrewer('qual', 'Dark2', 5));

load('params')
load('low_params')
params.d = 0; 

load('param_estimates_full');

factor = (param_estimates(1) + param_estimates(4))/param_estimates(1); 
params.cv =param_estimates(5); 

load('param_estimates_low.mat');
params.q_PI =  param_estimates(1);  
params.qT = param_estimates(2); 
params.beta_PI = param_estimates(3); 
params.beta_SI = factor*param_estimates(3); 

params.duration_cut_off = log10(1500);  

parameters = [params.time_start, params.time_end, params.Xinit, params.Yinit, params.Vinit...
params.Ninit, params.alpha, params.omega, params.d, 1e-10, 20, 1e-3, params.deltaT,params.dT,params.Tinit,1e-3, 5.5, 0, 0.5];

x = 1; y = 2;

%  %primary infection
params.time_end = 17; 
params.q = params.q_PI;
params.beta = params.beta_PI;

[T_PI, Y_PI] = ode45(@(t,y)PI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);  
V = log10(Y_PI(:,3));

params.beta = params.beta_PI - 2*params.cv*params.beta_PI;

[~, Y2_PI] = ode45(@(t,y)PI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit ]);  

params.beta = params.beta_PI + 2*params.cv*params.beta_PI; 
[~, Y3_PI] = ode45(@(t,y)PI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit ]);  

params.time_end = 15; 
params.beta = params.beta_SI;  
params.q = 0;  
params.qT_SI = params.qT; 

[T_SI, Y_SI] = ode45(@(t,y)SI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);  
V = log10(Y_SI(:,3));

params.beta = params.beta_SI - 2*params.cv*params.beta_SI; 
[~, Y2_SI] = ode45(@(t,y)SI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);  

params.beta = params.beta_SI + 2*params.cv*params.beta_SI; 
[~, Y3_SI] = ode45(@(t,y)SI(t, y, params),params.time_start:.01:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit params.Tinit]);  


%finding IP_js
params.IP = 5.9;
params.sigma = 0.2;

p =[0.001, 0.999];
crit = logninv(p,log(params.IP),params.sigma);
LB = round(crit(1), 2); 
UB = round(crit(2), 2); 
increment = (crit(2) - crit(1))/100; 
IP = LB:increment:UB; 

%PI
parameters(10) = params.beta_PI; %beta
parameters(11) = params.kappa; %kappa
parameters(12) = params.q_PI; 
parameters(13) = 0; %deltaT
parameters(16) = 0; %qT
parameters(5) = params.Vinit;

[T, Y] = eulers_method(0.001, parameters);

check = min(Y(3,:));
   
V = zeros(1, size(Y, 2));  
Ttemp = V; 

if check > 0 
     V = log10(Y(3, :));
     a = 0.001;
    Ttemp = round((1/a).*T);
end


data = params.data_PI_low; 
time = params.time_PI_low; 
temp_LOD = params.LOD_PI_low; 

subplot(x,y,1)
for j = 1:size(data, 1) %for each individual           
     
    index = isnan(data(j,:));
    time_temp = time(j,:);
    data_temp  = data(j,:);
    data_temp(isnan(data_temp)) = [];
    time_temp(index) = [];
    LOD =  temp_LOD(j,1);  
    
    counter = 0;

    for i = 1:length(IP)
        PI_posterior(j,i) =  find_IP_indiv(parameters, time_temp, data_temp, LOD, V, Ttemp, IP(i)) + log(lognpdf(IP(i),log(params.IP), params.sigma));
        
        counter = counter +  exp(PI_posterior(j,i));
        cdf_PI(j,i) = counter;
    end
   
    
   temp = sum(exp(PI_posterior(j,:))); 
   cdf_PI(j,:) = cdf_PI(j,:)./temp;
end

for j = 1:size(data, 1)
    val = unifrnd(0,1);
    vec = find(cdf_PI(j,:) <= val); 
    IP_val_PI(j) = IP(vec(end)); 
end

for j = 1:size(data, 1) 
   plot(params.time_PI_low(j,:) + IP_val_PI(j), params.data_PI_low(j,:),'k-.', 'LineWidth',1);
      hold on; 
end
     plot(0:15, log10(1500).*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','--'); 
    hold on;
    plot(0:15, log10(15000).*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','--');
    hold on; 
    plot(T_PI, log10(Y_PI(:,3)), 'Color', colors(3,:), 'LineWidth', 2);
    hold on;
    plot(T_PI, log10(Y2_PI(:,3)), 'Color', colors(3,:), 'LineWidth', 2, 'LineStyle', '--');
    hold on; 
    plot(T_PI, log10(Y3_PI(:,3)), 'Color', colors(3,:), 'LineWidth', 2, 'LineStyle', ':');
text(0.8,0.95,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 16)
ylabel('viral load (log genome copies/ml)'); title('1^° infection'); ylim([2.5,9]);
xlabel('time since infection (days)')
set(gca, 'FontSize', 16)
xlim([0, 15]);

%secondary infections
p = [0.001, 0.999];
crit = logninv(p, log(params.IP), params.sigma);
LB = round(crit(1), 2); 
UB = round(crit(2), 2); 
IP = LB:.01:UB; 

parameters(15) = 1e5; 
parameters(10) = params.beta_SI; %beta
parameters(11) = 5; %kappa
parameters(12) = 0; %q

parameters(13) = 1e-6; %deltaT
parameters(16) = params.qT_SI; %qT

[T, Y] = eulers_method(0.001, parameters);

check = min(Y(3,:));
   
V = zeros(1, size(Y, 2));  
Ttemp = V; 

if check > 0 
     V = log10(Y(3, :));
     a = 0.001;
    Ttemp = round((1/a).*T);
end

data = params.data_SI_low; 
time = params.time_SI_low; 
temp_LOD = params.LOD_SI_low; 

subplot(x,y,2)
for j = 1:size(data, 1) %for each individual           
     
    index = isnan(data(j,:));
    time_temp = time(j,:);
    data_temp  = data(j,:);
    data_temp(isnan(data_temp)) = [];
    time_temp(index) = [];
    LOD =  temp_LOD(j,1);  
   
    counter = 0;

    for i = 1:length(IP)
        SI_posterior(j,i) =  log(lognpdf(IP(i),log(params.IP), params.sigma)) + ... 
        find_IP_indiv(parameters, time_temp, data_temp, LOD, V, Ttemp, IP(i));

        counter = counter +  exp(SI_posterior(j,i));
        cdf_SI(j,i) = counter;
    end
    
    temp = sum(exp(SI_posterior(j,:))); 
    cdf_SI(j,:) = cdf_SI(j,:)./temp;
end

%sampling
for j = 1:size(data, 1)
    val = unifrnd(0,1);
    vec = find(cdf_SI(j,:) <= val); 
    IP_val_SI(j) = IP(vec(end)); 
end
 
for j = 1:size(data, 1)
   plot(params.time_SI_low(j,:) + IP_val_SI(j), params.data_SI_low(j,:),'k-.', 'LineWidth',1);
   hold on; 
end
plot(0:15, log10(1500).*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','--'); 
hold on;
plot(0:15, log10(15000).*ones(16,1), 'color', [0.5 0.5 0.5], 'LineStyle','--'); 
hold on; 
hold on; 
    plot(T_SI, log10(Y_SI(:,3)), 'Color', colors(4,:), 'LineWidth', 2);
    hold on;
    plot(T_SI, log10(Y2_SI(:,3)), 'Color', colors(4,:), 'LineWidth', 2, 'LineStyle', '--');
    hold on; 
    plot(T_SI, log10(Y3_SI(:,3)), 'Color', colors(4,:), 'LineWidth', 2, 'LineStyle', ':');
text(0.8,0.95,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 16)
title('2^° infection');  ylim([2.5, 9]); xlim([0, 15]);
xlabel('time since infection (days)') 
set(gca, 'FontSize', 16)

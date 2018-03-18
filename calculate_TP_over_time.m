 function TP = calculate_TP_over_time(viral_production_rate, s, params, infection_type)

%s = serotype (1-4)

if s == 1
    [final_params, ~] = fit_DENV_viral_load(1);  
elseif s == 2
    [final_params, ~] = fit_DENV_viral_load(2);  
elseif s ==3 
    [final_params, ~] = fit_DENV_viral_load(3);  
elseif s == 4
    [final_params, ~] = fit_DENV_viral_load(4);
end

mu = final_params(1); 
sigma = final_params(2); 
dt = 0.05;

params.omega = viral_production_rate;

if infection_type == 1
    [~, Y] = ode45(@(t,y)PI(t, y, params),params.time_start:dt:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);   
else
    [~, Y] = ode45(@(t,y)SI(t, y, params),params.time_start:dt:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit, params.Tinit]);   
end
V = log10(Y(:,3)); 
TP = V;
TP(1) = 1./(1 + exp(-sigma.*V(1) - mu));

for i = 2:length(V)
    TP(i) =  TP(i-1) + 1./(1 + exp(-sigma.*V(i) - mu))*dt;
end

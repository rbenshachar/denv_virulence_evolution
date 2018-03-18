function [TP, max_V, max_clearance] = calculate_SI_fitness_disease_alt(parameter_vector, s, params, dt)

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

p_hm = 1; m= 1; b = 1; 

mu = final_params(1); 
sigma = final_params(2);

TP = zeros(length(parameter_vector), 1);
duration = TP; 
max_V = TP; 
max_clearance = TP; 

for i = 1:length(parameter_vector)
 
   params.omega = parameter_vector(i);
   if i > 20 
       params.time_end = 50; 
   end

   [T, Y] = ode45(@(t,y)SI_alt(t, y, params),params.time_start:dt:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit, params.Tinit]);   
   V = log10(Y(:,3));  
   
   temp = find(V > params.duration_cut_off);
   
   if isempty(temp)
       temp = 2; 
   end
   
   V = V(temp(1):temp(end)); 
   T = T(temp(1):temp(end));
    
   q =  1./(1 + exp(-sigma.*V - mu));

   num_integral = 0; 
   for j = 1:length(T)
       num_integral = num_integral + q(j).*dt; 
   end

   duration(i) = T(end) - T(1); 
   
   max_V(i) = max(V);

   TP(i) = num_integral.*b^2.*m.*p_hm;

    temp = 0; 
    
    new_temp = find(V == max(V));  
    
    for z = new_temp(1):(length(V) - 1)
        new = (V(z) - V(z+1))/0.5;
        if new > temp
            temp = new; 
        end
 
    end  

    max_clearance = temp;
    
end

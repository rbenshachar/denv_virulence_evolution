 function [max_V, duration, TP, V, max_clearance] = calculate_PI_fitness_disease(parameter_vector, s, params, dt)

%s = serotype (1-4)
b=1; m=1; 
TP = zeros(length(parameter_vector), 1);
max_clearance = TP; 

max_V = TP; 
duration = TP;  

if s == 1
    final_params = fit_DENV_viral_load(1);  
elseif s == 2
    final_params = fit_DENV_viral_load(2); 
elseif s ==3 
    final_params = fit_DENV_viral_load(3); 
elseif s == 4
    final_params = fit_DENV_viral_load(4);
end

p_hm = 1; 
mu = final_params(1); 
sigma = final_params(2); 

detection = round(rand(1)*24); 

if detection == 0
    detection = 1; 
end

for i = 1:length(parameter_vector)
    if i > 50 
        params.time_end = 50; 
    end
   
    params.omega = parameter_vector(i);
 
   [T, Y] = ode45(@(t,y)PI(t, y, params),params.time_start:dt:params.time_end, [params.Xinit  params.Yinit params.Vinit params.Ninit]);   
   V = log10(Y(:,3)); 
    
   temp = find(V > params.duration_cut_off);
   
   if isempty(temp)
       temp = 2; 
   end
   
   V = V(temp(1):temp(end)); 
   T = T(temp(1):temp(end));
   
    noise_PI = zeros(length(V),1);
    l = 1;   

       for j = 1:length(noise_PI)  
          if  j == (detection + 12*(l-1))
             noise_PI(j) = normrnd(0,params.sigmae); 
             new_vals(l) = V(j) + noise_PI(j); 
             l = l + 1;  
           elseif j == 1
               noise_PI(1) = 0; 
           else
               noise_PI(j) = noise_PI(j-1);
           end 
       end
       
        if temp == 2
            new_vals = 1; 
        end

         V = V + noise_PI; 

       q =  1./(1 + exp(-sigma.*V - mu));

       num_integral = 0; 
       for j = 1:length(T)
           num_integral = num_integral + q(j).*dt; 
       end 

       max_V(i) = max(V); 

       TP(i) = num_integral.*b^2.*m.*p_hm;

       duration(i) = T(end) - T(1); 

        temp = 0; 
        new_temp = find(new_vals == max(new_vals)); 
        
        for z = new_temp(1):(length(new_vals) - 1)
            new = (new_vals(z) - new_vals(z+1))/0.5;
            if new > temp
                temp = new; 
            end
        end

        max_clearance = temp;
        
     if new_vals == 1
        TP(i) = 0; 
        max_V(i) = 0; 
    end
TP(i)
end
   

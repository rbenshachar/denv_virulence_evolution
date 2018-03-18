function value = new_likelihood_simpson_DE(parameters, time_temp, data_temp, LOD, V, Ttemp, x, f, c)

%taking log of likelihood
%increment is NOT on log-scale

increment = x(3) - x(1);
temp = zeros(1, length(x));

simulation_current_new = zeros(length(x), length(data_temp));  

parameters(17) = x(1); 

[simulation_current_new(1,:), T_vals] = model_fit(parameters, time_temp, V, Ttemp);
   
initial_IP = x(1); 

for i = 2:length(x)
    delta = x(i) - initial_IP; 
    T_new = T_vals + round(delta/.001);
    simulation_current_new(i,:) = V(T_new);
end

data_temp_new = repmat(data_temp, (length(x)), 1);
c1 = repmat(c, (length(x)), 1);
c2 = ones(size(c1)) - c1; 

temp(1:end) = get_likelihood(data_temp_new, LOD, simulation_current_new, parameters(19), c1, c2); 

left_f = f(1:2:length(f) - 1); 
left_temp = temp(1:2:length(temp) - 1); 

midpt_f = f(2:2:length(f)); 
midpt_temp = temp(2:2:length(temp)); 

right_f = f(3:2:length(f)); 
right_temp = temp(3:2:length(temp));  

simpson = (1/6).*increment.*(left_f.*left_temp + 4.*midpt_f.*midpt_temp + right_f.*right_temp);

value = sum(simpson);
   

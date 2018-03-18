function value = find_IP_indiv(parameters, time_temp, data_temp, LOD, V, Ttemp, x)

      c = zeros(1, length(data_temp)); 
   for i = 1:length(data_temp) %individual data points       
        if data_temp(i) <= LOD      
             c(i) = 0;  
        else
             c(i) = 1; 
        end
   end
  
  parameters(17) = x;
  parameters(19) = 0.5; 
 
  simulation_current = model_fit(parameters, time_temp, V, Ttemp);
 
  
  temp = get_likelihood(data_temp, LOD, simulation_current, parameters(19), c, 1-c); %likelihood not on log scale
 
 value = log(temp);

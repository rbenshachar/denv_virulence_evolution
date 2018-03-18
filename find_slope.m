function value = find_slope(time_temp, data_temp, LOD, x)

   %function to find slope of viral clearance for an individual
   index = isnan(data_temp);
   data_temp(isnan(data_temp)) = [];
   time_temp(index) = [];


   c = zeros(1, length(data_temp)); 
   for i = 1:length(data_temp) %individual data points       
        if data_temp(i) <= LOD      
             c(i) = 0;  
        else
             c(i) = 1; 
        end
   end
 
  Ttemp = 0:.001:25;
  Ttemp = round(Ttemp, 3); 
  model = x(1) + Ttemp.*x(2);
 
  [simulation_current, ~] = model_fit_no_IP(time_temp, model, Ttemp);
  temp = get_likelihood(data_temp, LOD, simulation_current, 0.5, c, 1-c); %likelihood not on log scale 
 
 value = -1*log(temp);


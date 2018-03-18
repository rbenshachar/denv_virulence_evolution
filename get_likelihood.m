function lik_value = get_likelihood(data_temp, LOD, simulation_temp, sigma_e, c1, c2)

%log-likelihood for each individual

   temp = (normcdf((LOD - simulation_temp)./sigma_e)).*(c2) + (normpdf((data_temp - simulation_temp)./sigma_e)).*(c1); 
   
   log_value = sum(log(temp), 2);
   lik_value = exp(log_value);
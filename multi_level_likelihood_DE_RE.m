function log_L = multi_level_likelihood_DE_RE(parameters, x, step_size, indiv, n)

lik_ind =  zeros(1,n); 

parameters(10) = x(1); %beta 
parameters(11) = x(2); %kappa
parameters(12) = x(3); %q 
IP_g = x(4);
parameters(13) = x(5); %deltaT
parameters(16) = x(6); %qT
parameters(18) = x(7); %sigma
parameters(5) = x(8); %Vinit 
parameters(19) = x(10); %sigmae
parameters(9) = x(11); %d

p =[0.025, 0.975];
crit = logninv(p,log(IP_g),parameters(18));

%simpson's rule for IPg
increment = (crit(2) - crit(1))/50; 
IP_range = crit(1):(increment/2):crit(2); 
f = lognpdf(IP_range, log(IP_g), parameters(18)); 

p2 = [0.025, 0.975]; %[0.25, 0.75];
crit_beta = norminv(p2,parameters(10),x(9));

%riemann sum for beta
increment_beta = (crit_beta(2) - crit_beta(1))/(step_size);
beta_range = crit_beta(1):(increment_beta/2):crit_beta(2);
 
f_beta = normpdf(beta_range, parameters(10), x(9)); 

simpson = zeros(1, length(f_beta)); 

V = zeros(length(f_beta), (parameters(2)-parameters(1))/0.001);  
Ttemp = zeros(1, (parameters(2)-parameters(1))/0.001); 

for k = 1:length(f_beta)            
        parameters(10) = beta_range(k);
        [T, Y] = eulers_method(0.001, parameters);

        check = min(Y(3,:));

        if check > 0 
             V(k,:) = log10(Y(3, :));
             a = 0.001;
             if k == 1
                  Ttemp = round((1/a).*T);
             end
        end  
end 

for j = 1:n %for each individual         
    data_temp = indiv(j).data;
    time_temp = indiv(j).time; 
    LOD =  indiv(j).LOD;
    
    for k = 1:length(beta_range)
        %integrating over all possible beta values using simpsons rule 
        temp = new_likelihood_simpson_DE(parameters, time_temp, data_temp, LOD, V(k,:), Ttemp, IP_range, f, indiv(j).c); 
        simpson(k) = temp;   
    end

   left_f = f_beta(1:2:length(f_beta) - 1);
   left_temp = simpson(1:2:length(simpson) - 1);
   
   midpt_f = f_beta(2:2:length(f_beta)); 
   midpt_temp = simpson(2:2:length(simpson)); 

   right_f = f_beta(3:2:length(f_beta)); 
   right_temp = simpson(3:2:length(simpson));  

   final_integral = (1/6).*increment_beta.*(left_f.*left_temp + 4.*midpt_f.*midpt_temp + right_f.*right_temp);
   lik_ind(j) = log(sum(final_integral));
   
   if isinf(lik_ind(j))
       lik_ind(j) = 0; 
   end
      
end

log_L = sum((lik_ind));

end
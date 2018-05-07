function dydt = simulate_dengue_ode_age(t, y, params_init, S_info, I_info)

[S_array, T_array, I_array, cumI_array] = UnVectorizeData(y, params_init);

dS_array = zeros(size(S_array)); 
dT_array = zeros(size(T_array));
dI_array = zeros(size(I_array)); 
dcumI_array = zeros(size(cumI_array));

% figure out rate of change:
I_array_ageless = sum(I_array, 3);

% first, compute force of infection for each serotype:
for i = 1:2 % only two serotypes co-circulating
    lambdaVals(i,1) = 0;
    for j = 1:(2^(params_init.n-1))
        if length(I_info(i,j).strain_history) == 0            
            lambdaVals(i,1) = lambdaVals(i,1) + params_init.beta_consti_matrix(i,1)*I_info(i,j).transFactor*(I_array_ageless(i,j)); % FOI now serotype-specific!!, also beta depends on PI, SI, postSI
        elseif length(I_info(i,j).strain_history) == 1
            lambdaVals(i,1) = lambdaVals(i,1) + params_init.beta_consti_matrix(i,2)*I_info(i,j).transFactor*(I_array_ageless(i,j)); 
        else % postsecondary infection
            lambdaVals(i,1) = lambdaVals(i,1) + params_init.beta_consti_matrix(i,3)*I_info(i,j).transFactor*(I_array_ageless(i,j)); 
        end
    end
    lambdaVals(i,1) = lambdaVals(i,1) + params_init.immigration;
end
lambdaVals(3:4,1) = 0;

% let's do S first: 
% S increases with births (NOT INCLUDED HERE - this is done at the daily level outside this function) 
% S decreases with new infections,
% S increases with T -> S

for i = 1:(2^params_init.n)
    
    % S decreases with new infections
    for j = 1:params_init.n
        if ~ismember(j, S_info(i).strain_history) % then susceptible to this strain.
            dS_array(:,i) = dS_array(:,i) - lambdaVals(j,1)*S_info(i).suscFactor*S_array(:,i)/params_init.N;  
        end
    end
    
    % S increases with T -> S
    dS_array(:,i) = dS_array(:,i) + params_init.delta*T_array(:,i);
    
end

% now let's do T: 
% T increases with I -> T
% T decreases with T -> S
% T increases or decreases with T -> T when there is clinical cross-protection

for i = 1:(2^params_init.n)
    
    % T increases with I -> T    
    for k = 1:params_init.n
        for l = 1:(2^(params_init.n - 1))
            new_infection_history = union(I_info(k,l).strain_history, k);
            if isequal(new_infection_history, S_info(i).strain_history)
                dT_array(:,i) = dT_array(:,i) + squeeze(params_init.v*I_array(k,l,:));
            end
        end
    end
        
    % T decreases with T -> S

    dT_array(:,i) = dT_array(:,i) - params_init.delta*T_array(:,i);
        
    if params_init.xImm_type == 1 % clinical cross-protection - T decreases
        for j = 1:params_init.n
            
            if ~ismember(j, S_info(i).strain_history) % then susceptible to this strain
                dT_array(:,i) = dT_array(:,i) - lambdaVals(j,1)*T_array(:,i)/params_init.N;
                
                % T increases with clinical cross-protection
                for k = 1:(2^params_init.n)
                    if isequal(union(S_info(i).strain_history, j), S_info(k).strain_history)
                        dT_array(:,k) = dT_array(:,k) + lambdaVals(j,1)*T_array(:,i)/params_init.N;
                        break;
                    end
                end
            end
        end
    end
end
    
% now let's do I: 
% I increases with S -> I
% I decreases with I -> T

for i = 1:params_init.n
    for j = 1:(2^(params_init.n - 1))
        
        % I increases with S -> I
        for k = 1:(2^params_init.n)
            if isequal(I_info(i,j).strain_history, S_info(k).strain_history)
                new_infections = lambdaVals(i,1)*S_info(k).suscFactor*S_array(:,k)/params_init.N;
                dI_array(i,j,:) = dI_array(i,j,:) + reshape(new_infections, 1, 1, params_init.n_age_classes);
            end 
        end
        
        dI_array(i,j,:) = dI_array(i,j,:) - params_init.v*I_array(i,j,:);
    
    end
end

% count cumulative infections upon recovery
dcumI_array = params_init.v*I_array; 

dydt = VectorizeData(dS_array, dT_array, dI_array, dcumI_array, params_init);

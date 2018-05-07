function dydt = simulate_dengue_ode_age_invasion(t, y, params_init, S_info, I_info)

[S_array, T_array, I_array, cumI_array] = UnVectorizeData(y, params_init);

dS_array = zeros(size(S_array)); 
dT_array = zeros(size(T_array));
dI_array = zeros(size(I_array)); 
dcumI_array = zeros(size(cumI_array));

% figure out rate of change:
I_array_ageless = sum(I_array, 3);

% first, compute force of infection for each serotype:
for i = 1:2  % only two serotypes
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

% dS_array stays at 0 -- invade when rare criteria
% dT_array stays at 0 -- invade when rare criteria

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

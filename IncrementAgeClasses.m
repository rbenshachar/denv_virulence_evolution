function [S_array, T_array, I_array] = IncrementAgeClasses(S_array, T_array, I_array, params)
   
% update arrays for mortalities:
for i = 1:params.n_age_classes
    S_array(i,:) = exp(-params.mu_rates(i))*S_array(i,:);
    if params.xImm_type ~= -1
        T_array(i,:) = exp(-params.mu_rates(i))*T_array(i,:);
    end
    I_array(:,:,i) = exp(-params.mu_rates(i))*I_array(:,:,i);
end

% now shift age classes: 
S_array(params.n_age_classes - 1,:) = S_array(params.n_age_classes - 1,:) + S_array(params.n_age_classes,:);
S_array(params.n_age_classes,:) = [];
S_array = [zeros(1, 2^params.n); S_array];

if params.xImm_type ~= -1
    T_array(params.n_age_classes - 1,:) = T_array(params.n_age_classes - 1,:) + T_array(params.n_age_classes,:);
    T_array(params.n_age_classes,:) = [];
    T_array = [zeros(1, 2^params.n); T_array];
end

I_array(:,:,params.n_age_classes - 1) = I_array(:,:,params.n_age_classes - 1) + I_array(:,:,params.n_age_classes);
I_array(:,:,params.n_age_classes) = [];

insert_I = zeros(params.n, 2^(params.n - 1), 1);
insert_I_vector = reshape(insert_I, params.n*(2^(params.n - 1)), 1);
I_array_vector = reshape(I_array, params.n*(2^(params.n - 1))*(params.n_age_classes-1), 1);
I_vector_new = [insert_I_vector; I_array_vector];
I_array_new = reshape(I_vector_new, params.n, 2^(params.n - 1), params.n_age_classes);
I_array = I_array_new;

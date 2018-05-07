function [S_array, T_array, I_array, cumI_array] = UnVectorizeData(y, params)

y_S = y(1:params.n_age_classes*(2^params.n));
S_array = reshape(y_S, params.n_age_classes, 2^params.n);
y(1:params.n_age_classes*(2^params.n)) = [];

if params.xImm_type == -1
    T_array = [];
else
    y_T = y(1:params.n_age_classes*(2^params.n));
    T_array = reshape(y_T, params.n_age_classes, 2^params.n);
    y(1:params.n_age_classes*(2^params.n)) = [];
end

y_I = y(1:params.n_age_classes*params.n*(2^(params.n-1)));
I_array = reshape(y_I, params.n, 2^(params.n-1), params.n_age_classes);
y(1:params.n_age_classes*params.n*(2^(params.n-1))) = [];

y_cumI = y(1:params.n_age_classes*params.n*(2^(params.n-1)));
cumI_array = reshape(y_cumI, params.n, 2^(params.n-1), params.n_age_classes);
y(1:params.n_age_classes*params.n*(2^(params.n-1))) = [];

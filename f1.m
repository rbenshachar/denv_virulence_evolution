function RMS = f1(n, v, x, parameters, indiv_PI, indiv_SI, fval, model)

if x > v(n)
    RMS = Inf; 
    return;
end

v(n) = x; 
if model == 1
    diff = MLE_likelihood_full(v, parameters, indiv_PI, indiv_SI) - fval;
else
    load('param_estimates_full.mat')
    a = (param_estimates(1) + param_estimates(4))/param_estimates(1);
    b = param_estimates(5);
    vals = [a,b];
    if model == 2
        diff = MLE_likelihood_subset_data(v, parameters, indiv_PI, indiv_SI, vals) - fval;
    elseif model == 3
        diff = MLE_likelihood_subset_data_innate(v, parameters, indiv_PI, indiv_SI, vals) - fval;
    elseif model == 4
        vals = [vals, 1e6];
        diff = MLE_likelihood_subset_data_diff_T(v, parameters, indiv_PI, indiv_SI, vals) - fval; 
    elseif model == 5    
        vals = [vals, 1e7];
        diff = MLE_likelihood_subset_data_diff_T(v, parameters, indiv_PI, indiv_SI, vals) - fval;
    end
end
x
RMS = abs(2 - diff)

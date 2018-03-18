function RMS = MLE_likelihood_subset_data_innate(v, parameters, indiv_PI, indiv_SI, vals)

close all; 

if min(v) < 0 
    RMS = Inf;
    return;
end 

x1 = cat(2, v(3), 5, v(1), 5.9,  0, 0, 0.2, 10^(-3),vals(2)*(v(3)), 0.5, 0); %PI
x2 = cat(2, vals(1)*v(3), 5, v(2), 5.9, 0, 0, 0.2, 10^(-3),vals(2)*(vals(1)*v(3)),0.5, 0); %SI

L_PI = multi_level_likelihood_DE_RE(parameters, x1, 50, indiv_PI, 11);
L_SI = multi_level_likelihood_DE_RE(parameters, x2, 50, indiv_SI, 90);

RMS = -L_PI - L_SI 
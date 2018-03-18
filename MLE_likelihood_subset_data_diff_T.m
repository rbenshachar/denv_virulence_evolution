function RMS = MLE_likelihood_subset_data_diff_T(v, parameters, indiv_PI, indiv_SI, vals)

close all; 

if min(v) < 0 
    RMS = Inf;
    return;
end 

x1 = cat(2, v(1), 5, v(2), 5.9,  0, 0, 0.2 ,10^(-3),vals(2)*(v(1)), 0.5, 1); %PI
x2 = cat(2, v(1)*vals(1), 5, 0,5.9, parameters(13), v(3), 0.2, 10^(-3),vals(2)*vals(1)*v(1),0.5, vals(3)); %SI

L_PI = multi_level_likelihood_DE_RE_2(parameters, x1, 50, indiv_PI, 11);
L_SI = multi_level_likelihood_DE_RE_2(parameters, x2, 50, indiv_SI, 90);
 
v
RMS =  -L_PI - L_SI 
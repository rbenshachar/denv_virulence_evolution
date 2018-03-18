function RMS = MLE_likelihood_full(v, parameters, indiv_PI, indiv_SI)

% close all; 
 
if min(v) < 0  
    RMS = Inf;
    return;
end 
%setting sigmae = 0.5
x1 = cat(2, v(1), 5, v(2), 5.9,  0, 0, 0.2,10^(-3),v(5)*(v(1)), 0.5, v(6)); %PI
x2 = cat(2, v(1) + v(4), 5, 0, 5.9, parameters(13), v(3),0.2, 10^(-3),v(5)*(v(1) + v(4)), 0.5, v(6)); %SI

L_PI = multi_level_likelihood_DE_RE(parameters, x1, 50, indiv_PI, 30);
L_SI = multi_level_likelihood_DE_RE(parameters, x2, 50, indiv_SI, 209);

%v
RMS = -L_PI - L_SI;
function final_params = fit_DENV_viral_load(k)

% logistic regression parameteres are from 
% Nguyen et al (2013) PNAS

if k == 1
    DV_OR = 3.93; DV_logV_at_p50percent = 6.5; % DV-1
elseif k == 2
     DV_OR = 5.48; DV_logV_at_p50percent = 6.3; % DV-2
elseif k == 3
    DV_OR = 2.61; DV_logV_at_p50percent = 7.7; % DV-3 
elseif k == 4
     DV_OR = 2.30; DV_logV_at_p50percent = 7.5; % DV-4
end

DV_beta1 = log(DV_OR); % since OR = exp(beta1)
DV_beta0 = -DV_beta1*DV_logV_at_p50percent;

final_params = [DV_beta0, DV_beta1];

end


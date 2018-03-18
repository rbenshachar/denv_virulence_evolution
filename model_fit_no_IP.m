function [model_V, T_vals] = model_fit_no_IP(time, model, Ttemp)

common_vals = ismember(Ttemp, time);

model_V =  model(common_vals);
T_vals = Ttemp(common_vals);

end


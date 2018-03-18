function [model_V, T_vals] = model_fit(parameters, time, V, Ttemp)

%RETURNS MODEL SIMULATED VIRAL LOAD LEVELS AT TIMEPOINTS WHERE HAVE DATA

%INPUTS:
%parameters
%time vector of data where viral load measurements were taken
%V: simulated viremia measurements
%Ttemp: simulated time vectors

%OUTPUTS:
%Simulated viral load levels at time points where have data
%corresponding time points

a = 0.001;  %# Define level of precision
timeRound = round((1/a).*(time + parameters(17))); % Round to appropriate precision

common_vals = ismember(Ttemp, timeRound); 

model_V =  V(common_vals);
 
T_vals = Ttemp(common_vals); 
 
end


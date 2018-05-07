function cumulativeI = RunDengueAgeSimulations_invasion(R0_PI_in, R0_SI_in, tStart, tEnd, dengue_infile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(dengue_infile);

params = params_init;

params.R0i_matrix = [R0_PI_in(1) R0_SI_in(1) 0; R0_PI_in(1) R0_SI_in(1) 0; 0 0 0; 0 0 0];
params.beta_consti_matrix = params.R0i_matrix*params.v;

[S_array, T_array, I_array, cumI_array] = UnVectorizeData(y_init(size(y_init,1),:)', params); cumI_array = zeros(size(cumI_array)); 
I_array_init = I_array;

t = []; y = [];
y_vals = VectorizeData(S_array, T_array, I_array, cumI_array, params);
[t, y] = ode45(@(t,y)simulate_dengue_ode_age_invasion(t, y, params, S_info_init, I_info_init), [tStart, tEnd], y_vals); 
[S_array, T_array, I_array, cumI_array] = UnVectorizeData(y(size(y, 1),:)', params);
 
cumulativeI = sum(sum(sum(cumI_array)));

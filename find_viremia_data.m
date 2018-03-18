function find_viremia_data()

close all;

[~, time_1_PI_DF, data_1_PI_DF, LOD_1_PI_DF] = viral_load_data(1,1,1);
[~, time_2_PI_DF, data_2_PI_DF, LOD_2_PI_DF] = viral_load_data(2,1,1);
[~, time_3_PI_DF, data_3_PI_DF, LOD_3_PI_DF] = viral_load_data(3,1,1);

[~, time_1_PI_DHF, data_1_PI_DHF, LOD_1_PI_DHF] = viral_load_data(1,1,2);
[~, time_2_PI_DHF, data_2_PI_DHF, LOD_2_PI_DHF] = viral_load_data(2,1,2);

[~, time_1_SI_DF, data_1_SI_DF, LOD_1_SI_DF] = viral_load_data(1,2,1);
[~, time_2_SI_DF, data_2_SI_DF, LOD_2_SI_DF] = viral_load_data(2,2,1);
[~, time_3_SI_DF, data_3_SI_DF, LOD_3_SI_DF] = viral_load_data(3,2,1);
[~, time_4_SI_DF, data_4_SI_DF, LOD_4_SI_DF] = viral_load_data(4,2,1);

[~, time_1_SI_DHF, data_1_SI_DHF, LOD_1_SI_DHF] = viral_load_data(1,2,2);
[~, time_2_SI_DHF, data_2_SI_DHF, LOD_2_SI_DHF] = viral_load_data(2,2,2);
[~, time_3_SI_DHF, data_3_SI_DHF, LOD_3_SI_DHF] = viral_load_data(3,2,2);
[~, time_4_SI_DHF, data_4_SI_DHF, LOD_4_SI_DHF] = viral_load_data(4,2,2);

params.time_PI_DF_1 = time_1_PI_DF; 
params.data_PI_DF_1 = data_1_PI_DF; 
params.LOD_PI_DF_1 = LOD_1_PI_DF; 

params.time_PI_DF_2 = time_2_PI_DF; 
params.data_PI_DF_2 = data_2_PI_DF; 
params.LOD_PI_DF_2 = LOD_2_PI_DF; 

params.time_PI_DF_3 = time_3_PI_DF; 
params.data_PI_DF_3 = data_3_PI_DF; 
params.LOD_PI_DF_3 = LOD_3_PI_DF; 

params.time_PI_DHF_1 = time_1_PI_DHF; 
params.data_PI_DHF_1 = data_1_PI_DHF; 
params.LOD_PI_DHF_1 = LOD_1_PI_DHF; 

params.time_PI_DHF_2 = time_2_PI_DHF; 
params.data_PI_DHF_2 = data_2_PI_DHF; 
params.LOD_PI_DHF_2 = LOD_2_PI_DHF; 

params.time_SI_DF_1 = time_1_SI_DF; 
params.data_SI_DF_1 = data_1_SI_DF; 
params.LOD_SI_DF_1 = LOD_1_SI_DF; 

%imputing missing LOD here
%LOD_2_SI_DF(21,1) = LOD_2_SI_DF(21,3);

% extra timept for one individual with no data
% taking this out
%time_2_SI_DF(21, 13) = NaN; 

params.time_SI_DF_2 = time_2_SI_DF; 
params.data_SI_DF_2 = data_2_SI_DF; 
params.LOD_SI_DF_2 = LOD_2_SI_DF; 

params.time_SI_DF_3 = time_3_SI_DF; 
params.data_SI_DF_3 = data_3_SI_DF; 
params.LOD_SI_DF_3 = LOD_3_SI_DF; 

params.time_SI_DF_4 = time_4_SI_DF; 
params.data_SI_DF_4 = data_4_SI_DF; 
params.LOD_SI_DF_4 = LOD_4_SI_DF; 

params.time_SI_DHF_1 = time_1_SI_DHF; 
params.data_SI_DHF_1 = data_1_SI_DHF; 
params.LOD_SI_DHF_1 = LOD_1_SI_DHF; 

params.time_SI_DHF_2 = time_2_SI_DHF; 
params.data_SI_DHF_2 = data_2_SI_DHF; 
params.LOD_SI_DHF_2 = LOD_2_SI_DHF; 

params.time_SI_DHF_3 = time_3_SI_DHF; 
params.data_SI_DHF_3 = data_3_SI_DHF; 
params.LOD_SI_DHF_3 = LOD_3_SI_DHF; 

params.time_SI_DHF_4 = time_4_SI_DHF; 
params.data_SI_DHF_4 = data_4_SI_DHF; 
params.LOD_SI_DHF_4 = LOD_4_SI_DHF; 

params.time_PI = cat(1, params.time_PI_DF_1, params.time_PI_DHF_1, params.time_PI_DF_2, ...
    params.time_PI_DHF_2, params.time_PI_DF_3); 
params.data_PI = cat(1, params.data_PI_DF_1, params.data_PI_DHF_1, params.data_PI_DF_2, ...
    params.data_PI_DHF_2, params.data_PI_DF_3); 
params.LOD_PI = cat(1, params.LOD_PI_DF_1, params.LOD_PI_DHF_1, params.LOD_PI_DF_2,...
    params.LOD_PI_DHF_2, params.LOD_PI_DF_3); 

params.time_SI = cat(1, params.time_SI_DF_1, params.time_SI_DHF_1, params.time_SI_DF_2,...
    params.time_SI_DHF_2, params.time_SI_DF_3, params.time_SI_DHF_3, ...
    params.time_SI_DF_4, params.time_SI_DHF_4); 
params.data_SI = cat(1, params.data_SI_DF_1, params.data_SI_DHF_1, params.data_SI_DF_2,...
    params.data_SI_DHF_2, params.data_SI_DF_3, params.data_SI_DHF_3, ...
    params.data_SI_DF_4, params.data_SI_DHF_4);
params.LOD_SI = cat(1, params.LOD_SI_DF_1, params.LOD_SI_DHF_1, params.LOD_SI_DF_2,...
    params.LOD_SI_DHF_2, params.LOD_SI_DF_3, params.LOD_SI_DHF_3,...
    params.LOD_SI_DF_4,  params.LOD_SI_DHF_4);

% params for simulations
params.time_start= 0;
params.time_end = 35;
params.kappa = 5; 
params.deltaT = 1.0e-06;
params.alpha = 1.0e-03;
params.omega = 1e4; 
params.Xinit = 10000000;
params.Yinit = 0;
params.Vinit = 1e-3; 
params.Ninit = 0; 
params.Tinit = 1e5; 
params.d = 0.07; 
params.dT = 0; 

save('params', 'params')
end
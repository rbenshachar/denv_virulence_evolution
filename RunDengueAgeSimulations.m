function void = RunDengueAgeSimulations(R0_PI, R0_SI, tStart, tEnd, dengue_outfile, dengue_infile, immunityType, durationImmunity)

load(dengue_infile);

params = params_init;

params.xImm_type = immunityType;
params.durationImmunity = durationImmunity;
params.delta = 1/params.durationImmunity; % delta is in yrs^-1

params.R0i_matrix = [R0_PI(1) R0_SI(1) 0; R0_PI(1) R0_SI(1) 0; 0 0 0; 0 0 0];
params.beta_consti_matrix = params.R0i_matrix*params.v;

%params.dt = 1/365; % 1 day increment
params.dt = 10/365; % 10 day increment
% the smaller the dt the better, but the slower

[S_array, T_array, I_array, cumI_array] = UnVectorizeData(y_init(size(y_init,1),:)', params); cumI_array = zeros(size(cumI_array)); 

t = []; y = [];

while 1
    
    % shift age classes up
    pre_N = sum(sum(S_array)) + sum(sum(T_array)) + sum(sum(sum(I_array)));
    [S_array, T_array, I_array] = IncrementAgeClasses(S_array, T_array, I_array, params);
    post_N = sum(sum(S_array)) + sum(sum(T_array)) + sum(sum(sum(I_array)));
    total_births_daily = (pre_N - post_N)*params.dt; % update births, by keeping population size constant over the timespan of a year - this is daily
    
    y_vals = VectorizeData(S_array, T_array, I_array, cumI_array, params);
    y_vals(1,1) = y_vals(1,1) + total_births_daily; % y_vals(1,1) corresponds to S_array(1,1)

    tFinal = tStart + 1;
    
    tStart_day = tStart
    while 1
        
        tFinal_day = tStart_day + params.dt; % simulate for 1 day
        [t_append, y_append] = ode45(@(t,y)simulate_dengue_ode_age(t, y, params, S_info_init, I_info_init), [tStart_day:(params.dt/2):tFinal_day], y_vals); 
        t = [t; t_append(2:length(t_append))]; y = [y; y_append(2:length(t_append),:)];
        tStart_day = tFinal_day;
        y_vals = y_append(length(t_append),:); y_vals(1,1) = y_vals(1,1) + total_births_daily;
        
        if (tFinal_day >= tFinal)
            break;
        end
    
    end   
    
    [S_array, T_array, I_array, cumI_array] = UnVectorizeData(y_append(size(y_append, 1),:)', params);
    tStart = tFinal;
        
    if (tFinal >= tEnd)
        params_init = params;
        t_init = t; y_init = y;
        save(dengue_outfile, 'params_init', 'S_info_init', 'I_info_init', 't_init', 'y_init');
        return;
    end
    
end

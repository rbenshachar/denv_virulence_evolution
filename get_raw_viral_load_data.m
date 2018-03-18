function [peakV, viral_clearance, peak] = get_raw_viral_load_data(k, s, i, method) 

%peakV is maximum viremia value
%viral clearance is viral clearance rate
%peak is if peak 0 or 1

%method = method for viral clearance
%method = 1: max viral clearance rate
%method = 2: inferred average sloped

%k = 1 (full dataset)/2 (subset of data)
%s = serotype 
%i = immune status
    
if k == 1 % full dataset
    [p1, time1, viral_load1, LOD1] = viral_load_data(s, i, 1); %DF
    
    t = 0; 
        
    if s == 1 || s == 2
        t = 1; 
        [p2, time2, viral_load2, LOD2] = viral_load_data(s, i, 2); %DHF
    end

    if i == 2 || t == 1 %have only primary DHF1 and primary DHF2
        if i == 2
            [p2, time2, viral_load2, LOD2] = viral_load_data(s, i, 2);
        end
        p = cat(2, p1, p2);
        time = cat(1,time1, time2); 
        viral_load = cat(1,viral_load1, viral_load2); 
        LOD = cat(1, LOD1, LOD2); 
	else
        p = p1; 
        time = time1; 
        viral_load = viral_load1; 
        LOD = LOD1; 
    end

else
    load('low_params') 
       
    if s == 1 && i == 1 %primary DENV-1
        time = params.time_PI_1_low;
        viral_load = params.data_PI_1_low;
    elseif s == 1 && i == 2
        time = params.time_SI_1_low; 
        viral_load = params.data_SI_1_low; 
    elseif s == 2 && i == 1
        time = params.time_PI_2_low; 
        viral_load = params.data_PI_2_low; 
    elseif s == 2 && i == 2
        time = params.time_SI_2_low; 
        viral_load = params.data_SI_2_low;
    elseif s == 3 && i == 1
        time = params.time_PI_3_low; 
        viral_load = params.data_PI_3_low;
    elseif s == 3 && i == 2
        time = params.time_SI_3_low; 
        viral_load = params.data_SI_3_low;
    elseif s == 4 && i == 1
        time = params.time_PI_4_low; 
        viral_load = params.data_PI_4_low;
    else
        time = params.time_SI_4_low; 
        viral_load = params.data_SI_4_low;
    end  
        
    p = time(:,1);  
end
    
peakV = zeros(1, length(p)); 
viral_clearance = peakV; 
peak = peakV; 
        
for j = 1:length(p)
            
    [peakV(j), index] = max(viral_load(j,:));
    
    if index(1) == 1
    	peak(j) = 0;
    else
        peak(j) = 1; 
    end
    
    if method == 1
        max_clearance = 0; 
        for l = 1:(length(time(j,:)) - 1)
            if isnan(time(j, l)) == 0 && isnan(time(j,l + 1)) == 0
                %rounding up to the nearest 0.5
                temp = (viral_load(j,l) - viral_load(j,l + 1))./(ceil(10*(time(j,l + 1) - time(j,l)))/10);
            end
            if temp(1) > max_clearance
                max_clearance = temp(1); 
            end
         end  

         viral_clearance(j) = max_clearance;
    else
        %Fitting slope of line
        lod = LOD(j,1);  

        data_temp =  viral_load(j,index:end);
        time_temp = round(time(j, index:end), 3);

         if size(data_temp) == 1 
             viral_clearance(j) = 0; 
         else
            x_init = [15, -2];  
            [param_estimates,  ~] = fminsearch(@(x)find_slope(time_temp, data_temp, lod, x), x_init);  
            viral_clearance(j) = -1*param_estimates(2); 
         end
    end          
end
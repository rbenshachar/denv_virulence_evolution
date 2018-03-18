function void = find_low_viremia_data(void)

close all; 

%30 individuals with PI 
%209 individuals with SI

%of lowest 50 %
%71 DENV-1, 21 DENV-2, 21 DENV-3, 6 DENV-4

load('params')

n1_PI = 18; n2_PI = 6; n3_PI = 6; n4_PI = 0; 
n1_SI = 124; n2_SI = 45; n3_SI = 33; n4_SI = 7; 

primary_peaks = zeros(size(params.data_PI, 1), 1); 
secondary_peaks = zeros(size(params.data_SI, 1), 1); 

%find 50% lowest PI
for i = 1:size(params.data_PI, 1)
    primary_peaks(i) = max(params.data_PI(i,:));
end

%peaks less than 8 logs - 37% of data
[PI_peaks, index_PI] = sort(primary_peaks); 
lowest_PI_peaks = index_PI(1:11);
% mean(PI_peaks(1:11))
% std(PI_peaks(1:11))
% median(PI_peaks(1:11))
% pause; 

temp1 = find(lowest_PI_peaks < (n1_PI + 1));
temp2 = find(lowest_PI_peaks > n1_PI & lowest_PI_peaks < (n1_PI + n2_PI + 1));
temp3 = find(lowest_PI_peaks > (n1_PI + n2_PI));

for i = 1:size(temp1)
    params.time_PI_1_low(i,:) = params.time_PI(lowest_PI_peaks(temp1(i)), :); 
    params.data_PI_1_low(i,:) = params.data_PI(lowest_PI_peaks(temp1(i)), :); 
    params.LOD_PI_1_low(i,:) = params.LOD_PI(lowest_PI_peaks(temp1(i)), :); 
end

for i = 1:size(temp2)
    params.time_PI_2_low(i,:) = params.time_PI(lowest_PI_peaks(temp2(i)), :); 
    params.data_PI_2_low(i,:) = params.data_PI(lowest_PI_peaks(temp2(i)), :); 
    params.LOD_PI_2_low(i,:) = params.LOD_PI(lowest_PI_peaks(temp2(i)), :);
end

for i = 1:size(temp3)
    params.time_PI_3_low(i,:) = params.time_PI(lowest_PI_peaks(temp3(i)), :); 
    params.data_PI_3_low(i,:) = params.data_PI(lowest_PI_peaks(temp3(i)), :); 
    params.LOD_PI_3_low(i,:) = params.LOD_PI(lowest_PI_peaks(temp3(i)), :);
end 

for i = 1:size(params.data_SI, 1)
    secondary_peaks(i) = max(params.data_SI(i,:));
end

%peaks < 8.5 magnitude, 43% of data
[peak_SI, index_SI] = sort(secondary_peaks); 
lowest_SI_peaks = index_SI(1:90); 

% mean(peak_SI(1:90))
% std(peak_SI(1:90))
% median(peak_SI(1:90))
%peak is used to check values

temp1 = find(lowest_SI_peaks < (n1_SI + 1));
temp2 = find(lowest_SI_peaks > n1_SI & lowest_SI_peaks < (n1_SI + n2_SI + 1));
temp3 = find(lowest_SI_peaks > (n1_SI + n2_SI) & lowest_SI_peaks < (n1_SI + n2_SI + n3_SI + 1));
temp4 = find(lowest_SI_peaks  > (n1_SI + n2_SI + n3_SI));

for i = 1:size(temp1)
    params.time_SI_1_low(i,:) = params.time_SI(lowest_SI_peaks(temp1(i)), :); 
    params.data_SI_1_low(i,:) = params.data_SI(lowest_SI_peaks(temp1(i)), :); 
    params.LOD_SI_1_low(i,:) = params.LOD_SI(lowest_SI_peaks(temp1(i)), :); 
end

for i = 1:size(temp2)
    params.time_SI_2_low(i,:) = params.time_SI(lowest_SI_peaks(temp2(i)), :); 
    params.data_SI_2_low(i,:) = params.data_SI(lowest_SI_peaks(temp2(i)), :); 
    params.LOD_SI_2_low(i,:) = params.LOD_SI(lowest_SI_peaks(temp2(i)), :);
end

for i = 1:size(temp3)
    params.time_SI_3_low(i,:) = params.time_SI(lowest_SI_peaks(temp3(i)), :); 
    params.data_SI_3_low(i,:) = params.data_SI(lowest_SI_peaks(temp3(i)), :); 
    params.LOD_SI_3_low(i,:) = params.LOD_SI(lowest_SI_peaks(temp3(i)), :);
end

for i = 1:size(temp4)
    params.time_SI_4_low(i,:) = params.time_SI(lowest_SI_peaks(temp4(i)), :); 
    params.data_SI_4_low(i,:) = params.data_SI(lowest_SI_peaks(temp4(i)), :); 
    params.LOD_SI_4_low(i,:) = params.LOD_SI(lowest_SI_peaks(temp4(i)), :);
end

params.time_PI_low = cat(1, params.time_PI_1_low, params.time_PI_2_low,...
    params.time_PI_3_low);
params.time_SI_low = cat(1, params.time_SI_1_low, params.time_SI_2_low,...
    params.time_SI_3_low, params.time_SI_4_low);
params.data_PI_low = cat(1, params.data_PI_1_low, params.data_PI_2_low,... 
    params.data_PI_3_low);
params.data_SI_low = cat(1, params.data_SI_1_low, params.data_SI_2_low,...
    params.data_SI_3_low, params.data_SI_4_low);
params.LOD_PI_low = cat(1, params.LOD_PI_1_low, params.LOD_PI_2_low,... 
    params.LOD_PI_3_low);
params.LOD_SI_low = cat(1, params.LOD_SI_1_low, params.LOD_SI_2_low,...
    params.LOD_SI_3_low, params.LOD_SI_4_low);

save('low_params', 'params')

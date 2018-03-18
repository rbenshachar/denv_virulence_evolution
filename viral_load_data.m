function [pat_num, time, viral_load, LOD] = viral_load_data(serotype, type_infection, CM)

% this file is formatted from the python file: format_data.py
M = xlsread('Clapham_data_formatted.xlsx');

j = 1; 

for i = 1:size(M, 1)
    if M(i,2) == serotype && M(i,4) == type_infection && M(i,5) == CM  && M(i,1)~= M(i-1,1)
        pat_num(j) = M(i,1);
        j = j+1; 
    end
end

time = NaN(length(pat_num), 16); 
viral_load = NaN(length(pat_num), 16);
LOD = NaN(length(pat_num), 16); 


for i = 1:length(pat_num)
    temp = find(M(:,1) == pat_num(i)); 
    time(i,1:length(temp)) = M(temp, 6); 
    viral_load(i,1:length(temp)) = log10(M(temp, 7));
    LOD(i,1:length(temp)) = log10(M(temp, 8)); 
    
    for j = 1:length(temp)
        if viral_load(i,j) < LOD(i,j) 
            viral_load(i,j) = LOD(i,j); 
        end
   end
end

if serotype == 2 && type_infection == 2 && CM == 1
    %imputing missing LOD here
    LOD(21,1) = LOD(21,3);
end
function IP_val = finding_IPs(parameters, data, time, temp_LOD, params, IP)

[T, Y] = eulers_method(0.001, parameters);

check = min(Y(3,:));
   
V = zeros(1, size(Y, 2));  
Ttemp = V; 

if check > 0 
     V = log10(Y(3, :));
     a = 0.001;
    Ttemp = round((1/a).*T);
end

for j = 1:size(data, 1) %for each individual           
     
    index = isnan(data(j,:));
    time_temp = time(j,:);
    data_temp  = data(j,:);
    data_temp(isnan(data_temp)) = [];
    time_temp(index) = [];
    LOD =  temp_LOD(j,1);  
    
    counter = 0;

    for i = 1:length(IP)
        PI_posterior(j,i) =  log(lognpdf(IP(i),log(5.9), params.sigma)) + ... 
        find_IP_indiv(parameters, time_temp, data_temp, LOD, V, Ttemp, IP(i));
        counter = counter +  exp(PI_posterior(j,i));
        cdf_PI(j,i) = counter;
    end
    
    temp = sum(exp(PI_posterior(j,:))); 
    cdf_PI(j,:) = cdf_PI(j,:)./temp;
end

for j = 1:size(data, 1)
    val = unifrnd(0,1);
    vec = find(cdf_PI(j,:) <= val); 
    IP_val(j) = IP(vec(end));
end
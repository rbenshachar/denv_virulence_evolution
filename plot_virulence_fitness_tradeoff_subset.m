function plot_virulence_fitness_tradeoff_subset(s)

%%%%%%%%%%%%%%%%%
%FIGURE 3

% s is serotype 
% plotting R0 as a function of peak V for primary and secondary

close all;

load('params');

colors = colormap(cbrewer('qual', 'Dark2', 5));

k = 100; 

load('param_estimates_full.mat')
cv = param_estimates(5);
f = (param_estimates(1) + param_estimates(4))/param_estimates(1); 

% %FIT TO PARTIAL DATASET
load('param_estimates_low.mat')
params.cv = cv; 
params.sigmae = 0; 
params.d = 0;  
params.duration_cut_off = log10(1500); 

params.q_PI = param_estimates(1); 
params.qT_SI = param_estimates(2); 
params.beta_PI = param_estimates(3); 
params.beta_SI = param_estimates(3)*f; 
params.sd_PI = params.beta_PI*params.cv; 
params.sd_SI = params.beta_SI*params.cv; 

phenotype = 1e4*(.01:.02:4); %1e4*(.01:.02:10);

TP_PI = zeros(k, length(phenotype)); 
disease_PI = zeros(k,length(phenotype));

TP_SI = zeros(k, length(phenotype)); 
disease_SI = zeros(k,length(phenotype));

% %PI
for i = 1:k
    %omega
    params.time_end = 150; 
    params.beta = normrnd(params.beta_PI, params.sd_PI);
    params.q = params.q_PI; 

    [disease_PI(i,:), ~, TP_PI(i,:), ~, ~] = calculate_PI_fitness_disease(phenotype, s, params, .01);
    i
end

  params.time_end = 200; 
    params.beta = params.beta_PI;
    params.q = params.q_PI;

    [disease_PI_mean, ~, TP_PI_mean, ~, ~] = calculate_PI_fitness_disease(phenotype, s, params, .01);

% %SI
for i = 1:k
    %omega
    params.time_end = 250; 
    params.beta = normrnd(params.beta_SI, params.sd_SI);
    params.Tinit = 1e5; 
    params.q = 0;  
    params.qT = params.qT_SI;

    [TP_SI(i,:), disease_SI(i,:), ~] = calculate_SI_fitness_disease(phenotype, s, params, .01);
    i
end

params.time_end = 250; 
params.beta = params.beta_SI;
params.Tinit = 1e5; 
params.q = 0;  
params.qT = params.qT_SI;

[TP_SI_mean, disease_SI_mean, ~] = calculate_SI_fitness_disease(phenotype, s, params, .01);

upper_TP_omega_PI = zeros(length(phenotype), 1)'; 
lower_TP_omega_PI = upper_TP_omega_PI; 
upper_disease_omega_PI = zeros(length(phenotype), 1)'; 
lower_disease_omega_PI = zeros(length(phenotype), 1)'; 

for j = 1:length(phenotype)
   upper_TP_omega_PI(j) = quantile(TP_PI(:,j), 0.975); 
   lower_TP_omega_PI(j) = quantile(TP_PI(:,j), 0.025); 
   upper_disease_omega_PI(j) =  quantile(disease_PI(:,j), 0.975); 
   lower_disease_omega_PI(j) = quantile(disease_PI(:,j), 0.025); 
end

upper_TP_omega_SI = zeros(length(phenotype), 1)'; 
lower_TP_omega_SI = upper_TP_omega_SI; 
upper_disease_omega_SI = zeros(length(phenotype), 1)'; 
lower_disease_omega_SI = zeros(length(phenotype), 1)'; 

for j = 1:length(phenotype)
   upper_TP_omega_SI(j) = quantile(TP_SI(:,j), 0.975); 
   lower_TP_omega_SI(j) = quantile(TP_SI(:,j), 0.025); 
   upper_disease_omega_SI(j) = quantile(disease_SI(:,j), 0.975); 
   lower_disease_omega_SI(j) = quantile(disease_SI(:,j), 0.025); 
end
 
% save('TP_PI_D1', 'TP_PI')
% save('TP_SI_D1', 'TP_SI')
% 
% save('disease_PI_D1', 'disease_PI')
% save('disease_SI_D1', 'disease_SI')

x = 3; y = 1;

subplot(x,y,1)
shadedplot(phenotype, lower_disease_omega_PI,  upper_disease_omega_PI, colors(1,:), 'none');
hold on; 
plot(phenotype, disease_PI_mean, 'k', 'LineWidth', 2);
hold on; 
shadedplot(phenotype, lower_disease_omega_SI, upper_disease_omega_SI, colors(2,:), 'none');
hold on; 
plot(phenotype, disease_SI_mean, 'k', 'LineWidth', 2); 
text(0.9,0.95,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
xlabel('viral production rate (genome copies/cell/day)'); ylabel('peak viral load (log genome copies/ml)')
set(gca, 'FontSize', 14)

subplot(x,y,2)
shadedplot(phenotype, lower_TP_omega_PI,  upper_TP_omega_PI, colors(1,:), 'none');
hold on; 
plot(phenotype, TP_PI_mean, 'k', 'LineWidth', 2);
hold on; 
shadedplot(phenotype, lower_TP_omega_SI, upper_TP_omega_SI, colors(2,:), 'none');
hold on; 
plot(phenotype, TP_SI_mean, 'k', 'LineWidth', 2); 
text(0.9,0.95,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
xlabel('viral production rate (genome copies/cell/day)'); ylabel('transmission potential')
set(gca, 'FontSize', 14)

%subplot c
z = 0.5;

y = disease_PI'; 
disease_PI_vec = y(:)';
disease_PI_vec = round(disease_PI_vec./z).*z; %round(disease_PI_vec/0.05)*0.05; 

y = R0_PI';
TP_PI_vec = y(:)';

[disease_PI_vec, index] = sort(disease_PI_vec); 

temp = find(disease_PI_vec == 0); 
if isempty(temp) == 0
    disease_PI_vec = disease_PI_vec((temp(end)+1):end); 
    index = index((temp(end) + 1):end); 
end

TP_PI_vec = TP_PI_vec(index); 
 
peakV_PI = round(((min(disease_PI_vec)):z:(max(disease_PI_vec)))./z).*z;
peakV_TP_PI = zeros(3, length(peakV_PI));   

for i = 1:length(peakV_PI)
    
    temp = TP_PI_vec(disease_PI_vec == peakV_PI(i));   
    peakV_TP_PI(1,i) = quantile(temp, 0.025); 
    peakV_TP_PI(2,i) = quantile(temp, 0.5); 
    peakV_TP_PI(3,i) = quantile(temp, 0.975); 
    
end

y = disease_SI'; 
disease_SI_vec = y(:)';
disease_SI_vec = round(disease_SI_vec./z).*z;

y = TP_SI';
TP_SI_vec = y(:)';

[disease_SI_vec, index] = sort(disease_SI_vec); 

temp = find(disease_SI_vec == 0); 
if isempty(temp) == 0
    disease_SI_vec = disease_SI_vec((temp(end)+1):end); 
    index = index((temp(end) + 1):end); 
end

TP_SI_vec = TP_SI_vec(index); 

peakV_SI = round(((min(disease_SI_vec)):z:(max(disease_SI_vec)))./z).*z;
peakV_TP_SI = zeros(3, length(peakV_SI));   

for i = 1:length(peakV_SI)
    
    temp = TP_SI_vec(disease_SI_vec == peakV_SI(i));   
    peakV_TP_SI(1,i) = quantile(temp, 0.025); 
    peakV_TP_SI(2,i) = quantile(temp, 0.5); 
    peakV_TP_SI(3,i) = quantile(temp, 0.975);
    
end

x = 3; y = 1;
subplot(x,y,3)
shadedplot(peakV_PI, peakV_TP_PI(1,:), peakV_TP_PI(3,:), colors(1,:), 'none');
hold on; 
plot(peakV_PI, peakV_TP_PI(2,:), 'k', 'LineWidth', 2); 
hold on; 
shadedplot(peakV_SI, peakV_TP_SI(1,:), peakV_TP_SI(3,:), colors(2,:), 'none');
hold on; 
plot(peakV_SI, peakV_TP_SI(2,:), 'k', 'LineWidth', 2); 
text(0.9,0.95,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
xlabel('peak viral load (log copies/ml)'); ylabel('transmission potential')
set(gca, 'FontSize', 14)



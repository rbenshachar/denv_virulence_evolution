function [p1, a1, p2, a2] = plot_viral_peak_clearance(method) 

%%FIGURE 1

%% method= 1: viral clearance is measured as max viral clearance (Fig. 1)
%% method= 2: fitting slope of viral clearance (additional method mentioned 
%% in Methods) 

close all; 
colors = colormap(cbrewer('qual', 'Paired', 8));
x = 1; y = 2; 

%get viral load data
find_viremia_data();

load('params')

%full dataset
%Primary infections
n = 1; 
for i = 1:3 %serotypes - no primary DENV-4      
    [peakV, slope,peak] = get_raw_viral_load_data(1, i, 1, method);   
    for j = 1:length(peak)
        if peak(j) == 1
            y1(n) = peakV(j); 
            x1(n) = slope(j);  
            n = n+ 1; 
        end
    end
end   

mdl = fitlm(y1', x1', 'linear'); 

% using cook's distance
outliers = find((mdl.Diagnostics.CooksDistance)>3*mean(mdl.Diagnostics.CooksDistance));
temp = 1:length(y1);
z = setdiff(temp, outliers);

y1 = y1(z); 
x1 = x1(z); 

subplot(x,y,1)
for j = 1:length(y1)
    plot(y1(j), x1(j) , 'o', 'Color', colors(1,:), 'MarkerSize', 10, 'MarkerFaceColor','none', 'LineWidth', 2);
    hold on;
end

f1 = fitlm(y1,x1, 'linear');
a1 = f1.Coefficients.Estimate; 
p1 = f1.Coefficients.pValue(2);
   
subplot(x,y,1)
v = 3:12;
plot(v, a1(2).*v + a1(1), 'k','LineWidth', 2);              
text(0.8,0.9,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
xlim([5, 12]); ylim([1, 8])
set(gca, 'FontSize', 14);
   
if method == 1
    xlabel('peak viral load (log genome copies/ml)'); ylabel('maximum daily viral clearance rate (log copies/ml/day)') 
else
    xlabel('peak viral load (log genome copies/ml)'); ylabel('average daily viral clearance rate (log copies/ml/day)') 
end

n = 1; 
for i = 1:4 %four serotypes
    [peakV, slope,peak] = get_raw_viral_load_data(1, i, 2, method);

    subplot(x,y,2)
    for j = 1:length(peak)
        if peak(j) == 1
            plot(peakV(j), slope(j), 'o', 'Color', colors(1,:), 'MarkerSize', 10, 'MarkerFaceColor','none', 'LineWidth', 2); 
            hold on; 
            y1(n) = peakV(j); 
            x1(n) = slope(j);
            n = n+ 1; 
        end
    end
end

mdl = fitlm(y1',x1', 'linear');

%using cook's distance
outliers = find((mdl.Diagnostics.CooksDistance)>3*mean(mdl.Diagnostics.CooksDistance));
temp = 1:length(y1);
z = setdiff(temp, outliers);

y1 = y1(z); 
x1 = x1(z); 

f2 = fitlm(y1',x1', 'linear');

a2 = f2.Coefficients.Estimate; 
p2 = f2.Coefficients.pValue(2);
   
subplot(x,y,2)
v = 3:12;
plot(v, a2(2).*v + a2(1), 'k','LineWidth', 2); 
if method  == 1 
    xlabel('peak viral load (log genome copies per ml)'); ylabel('maximum daily viral clearance rate (log genome copies per ml per day)') 
else
    xlabel('peak viral load (log genome copies per ml)'); ylabel('average daily viral clearance rate (log genome copies per ml per day)') 
end
set(gca, 'FontSize', 14);
text(0.8,0.9,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 14)
xlim([5, 12]); ylim([1, 8]) 

end

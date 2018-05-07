function void = PlotInfecteds(infile)

load(infile)

[ntimes, ncols] = size(y_init);

for t = 1:ntimes
    [S_array, T_array, I_array, cumI_array] = UnVectorizeData(y_init(t,:)', params_init); 
 
    S_all(t) = sum(sum(S_array));
    T_all(t) = sum(sum(T_array));
    I_age_less = sum(I_array, 3);
    n_primary_infections(t) = I_age_less(1,1) + I_age_less(2,1);
    n_secondary_infections(t) = I_age_less(1,2) + I_age_less(2,2);
end
subplot(2,1,1); plot(t_init, n_primary_infections, 'b'); hold on;
plot(t_init, n_secondary_infections, 'r');  xlabel('time (yrs)'); ylabel('# infecteds'); legend('primary infections', 'secondary infections');

subplot(2,1,2); plot(t_init, S_all, 'm'); hold on;
plot(t_init, T_all, 'k');  xlabel('time (yrs)'); ylabel('# infecteds'); legend('all S', 'all T');
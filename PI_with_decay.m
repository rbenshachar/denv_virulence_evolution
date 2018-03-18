function ydot = PI_with_decay(t, y_in, params)

x = y_in(1);
y = y_in(2);
v = y_in(3); 
n = y_in(4); 

dxdt = - params.beta*x*v;
dvdt = params.omega*y - params.kappa*v;
dndt = params.q*y - params.d*n;  
if t > 5 && v < 1500
    dydt = params.beta*x*v - params.alpha*n*y - params.Tinit*y*params.deltaT*exp(params.delta*(t - 4));
else
   dydt = params.beta*x*v - params.alpha*n*y; 
end

ydot = [dxdt dydt dvdt dndt]'; 
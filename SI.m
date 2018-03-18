function ydot = SI(t, y_in, params)

x = y_in(1);
y = y_in(2);
v = y_in(3); 
n = y_in(4); 
T = y_in(5);  

dxdt = - params.beta*x*v;
dydt = params.beta*x*v - params.alpha*n*y- params.deltaT*y*T; 
dvdt = params.omega*y - params.kappa*v;
dndt = params.q*y - params.d*n;  
dTdt = params.qT*y*T;

ydot = [dxdt dydt dvdt dndt dTdt]';
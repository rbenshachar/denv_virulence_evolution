function [T, Y] =  eulers_method_2(h, parameters)

    time_start = parameters(1); 
    time_end = parameters(2); 
    Xinit = parameters(3);  
    Yinit = parameters(4); 
    Vinit = parameters(5); 
    Ninit = parameters(6); 
    alpha  = parameters(7); 
    omega = parameters(8); 
    d = parameters(9); 
    beta = parameters(10); 
    kappa = parameters(11); 
    q = parameters(12);
    deltaT = parameters(13); 
    dT = parameters(14); 
    Tinit = parameters(15); 
    qT = parameters(16); 
    a = parameters(20); 
    
    n = (time_end - time_start)/h; 

    T = zeros(1, n); 
    Y = zeros(4, n); 
    T_cells = zeros(1, n); 

    Y(1,1) = Xinit; 
    Y(2,1) = Yinit; 
    Y(3,1) = Vinit; 
    Y(4,1) = Ninit; 
    T_cells(1,1) = Tinit; 

    for i = 2:n

        T(i) = T(i-1) + h; 

        Y(1,i) = Y(1, i-1) + h*(- beta*Y(1,i-1)*Y(3, i-1));
        Y(2,i) = Y(2, i-1) + h*(beta*Y(1,i-1)*Y(3, i-1) - alpha*Y(4, i-1)*Y(2, i-1) - deltaT*T_cells(1, i-1)*Y(2, i-1)); 
        Y(3,i) = Y(3, i-1) + h*(omega*Y(2, i-1) - kappa*Y(3, i-1)); 
        Y(4,i) = Y(4, i-1) + h*(q*Y(2, i - 1) - d*Y(4, i-1)); 
        
        T_cells(1,i) = T_cells(1, i-1) + h*(qT*T_cells(1, i-1)*Y(2, i-1)/(a + Y(2, i-1))); 
   
    end

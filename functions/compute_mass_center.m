function [x_cm] = compute_mass_center(material,rhoA,naca_2dt,c)
    
    x_cm = 0;
    % Leading edge
    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,0,0.0075,1000,c);
    integrand = airfoil_points.*x;
    x_cm = x_cm + material.Nylon.rho*2*trapz(x,integrand);%*1e-6;

    % Plate
    x = linspace(0.020,0.080,20);
    integrand = x*0.00225;
    x_cm = x_cm + material.Al.rho*trapz(x,integrand);%*1e-6;

    % Trailing edge
    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,0.095,0.100,1000,c);
    integrand = airfoil_points.*x;
    x_cm = x_cm + material.Nylon.rho*2*trapz(x,integrand);%*1e-6;
    
    x_cm = x_cm/rhoA;
end
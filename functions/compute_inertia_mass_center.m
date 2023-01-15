function [I_cm] = compute_inertia_mass_center(x_cm,material,naca_2dt,c)
    
    I_cm = 0;
    
    % Leading edge
    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,0,0.0075,1000,c);
    integrand = airfoil_points.*(x-x_cm).^2;
    I_cm = I_cm + material.Nylon.rho*2*trapz(x,integrand);%*1e-6;

    % Plate
    x = linspace(0.020,0.080,20);
    integrand = 0.00225*(x-x_cm).^2;
    I_cm = I_cm + material.Al.rho*trapz(x,integrand);%*1e-6;

    % Trailing edge
    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,0.095,0.100,1000,c);
    integrand = airfoil_points.*(x-x_cm).^2;
    I_cm = I_cm + material.Nylon.rho*2*trapz(x,integrand);%*1e-6;

end
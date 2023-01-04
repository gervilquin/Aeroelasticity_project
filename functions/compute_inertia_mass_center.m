function [I_cm] = compute_inertia_mass_center(x_cm,material,naca_2dt,c)
    
    I_cm = 0;
    % Leading edge
    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,0,7.5,1000,c);
    integrand = airfoil_points.*(x-x_cm).^2;
    I_cm = I_cm + material.Nylon.rho*2*trapz(x,integrand)*1e-6;

    % Plate
    x = linspace(20,80,20);
    integrand = 2.25*(x-x_cm).^2;
    I_cm = I_cm + material.Al.rho*trapz(x,integrand)*1e-6;

    % Trailing edge
    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,95,100,1000,c);
    integrand = airfoil_points.*(x-x_cm).^2;
    I_cm = I_cm + material.Nylon.rho*2*trapz(x,integrand)*1e-6;

end
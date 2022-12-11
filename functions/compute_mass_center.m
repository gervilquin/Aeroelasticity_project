function [x_cm] = compute_mass_center(material,rhoA,naca_2dt,c)
    x_cm = 0;
    % Leading edge
    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,0,7.5,20,c);
    integrand = airfoil_points.*x;
    x_cm = x_cm + material.Al.rho*2*trapz(x,integrand)*1e-6;

    % Plate
    x = linspace(20,80,20);
    integrand = x*2.25;
    x_cm = x_cm + material.Al.rho*trapz(x,integrand)*1e-6;

    % Trailing edge
    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,95,100,20,c);
    integrand = airfoil_points.*x;
    x_cm = x_cm + material.Al.rho*2*trapz(x,integrand)*1e-6;

    x_cm = x_cm/rhoA;
end
function [rhoA] = compute_avg_density(material,naca_2dt,c)
% ------------------------
    rhoA = 0;

    % Leading edge
    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,0,0.0075,20,c);
    rhoA = rhoA + material.Nylon.rho*2*trapz(x,airfoil_points);%*1e-6;
    
    % Plate
    rhoA = rhoA + material.Al.rho*(60*2.25)*1e-6;

    % Trailing edge
    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,0.095,0.100,20,c);
    rhoA = rhoA + material.Nylon.rho*2*trapz(x,airfoil_points);%*1e-6;
end
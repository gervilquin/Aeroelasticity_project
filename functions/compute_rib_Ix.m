function Ix_rib = compute_rib_Ix(rho_r,naca_2dt,c)

    slot_start = 35;
    slot_end = 65;
    slot_height = 5;

    Ix = 0;

    % Leading edge
    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,7.5,slot_start,1000,c);
    integrand = airfoil_points.^3;
    Ix = Ix + 2*trapz(x,integrand)*1e-6;

    % Section with the slot
    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,slot_start,slot_end,1000,c);
    slot_points = slot_height/2*ones(1,length(airfoil_points));
    integrand = airfoil_points.^3-slot_points.^3;
    Ix = Ix + 2*trapz(x,integrand)*1e-6;
        

    % Trailing edge
    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,slot_end,95,1000,c);
    integrand = airfoil_points.^3;
    Ix = Ix + 2*trapz(x,integrand)*1e-6;

    Ix_rib = Ix*rho_r;

end
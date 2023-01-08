function xcm_r = compute_rib_xcm(A_r,naca_2dt,c)
    
    slot_start = 35;
    slot_end = 65;
    slot_height = 5;

    x_cm = 0;

    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,7.5,slot_start,1000,c);
    integrand = airfoil_points.*x;
    x_cm = x_cm + 2*trapz(x,integrand)*1e-6;

    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,slot_start,slot_end,1000,c);
    integrand = (airfoil_points-ones(1,length(airfoil_points))*slot_height/2).*x;
    x_cm = x_cm + 2*trapz(x,integrand)*1e-6;

    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,slot_end,95,1000,c);
    integrand = airfoil_points.*x;
    x_cm = x_cm + 2*trapz(x,integrand)*1e-6;

    xcm_r = x_cm/A_r;
end
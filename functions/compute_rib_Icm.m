function Icm_r = compute_rib_Icm(rho_r,xcm_r,naca_2dt,c)

    slot_start = 35;
    slot_end = 65;
    slot_height = 5;

    I_cm = 0;

    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,7.5,slot_start,1000,c);
    integrand = airfoil_points.*(x-xcm_r).^2;
    I_cm = I_cm + 2*trapz(x,integrand)*1e-6;

    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,slot_start,slot_end,1000,c);
    integrand = (airfoil_points-ones(1,length(airfoil_points))*slot_height/2).*(x-xcm_r).^2;
    I_cm = I_cm + 2*trapz(x,integrand)*1e-6;

    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,slot_end,95,1000,c);
    integrand = airfoil_points.*(x-xcm_r).^2;
    I_cm = I_cm + 2*trapz(x,integrand)*1e-6;

    Icm_r = I_cm*rho_r;


end
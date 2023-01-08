function A_rib = compute_rib_area(naca_2dt,c)
    
    [airfoil_points,x] = compute_symmetric_airfoil(naca_2dt,7.5,95,100,c);
    A_rib = 2*trapz(x,airfoil_points)*1e-6;

    A_slot = 5*30*1e-6;
    
    A_rib = A_rib - A_slot;
end
function Icm_r = ComputeInertiaCMRib(rho_r,xcm_r,t,c)
    I_cm = 0;

    % Computation without slot
    [x,airfoil_points] = CoordinatesSymmetricAirfoil(t,0.0075,0.095,200,c,'TE');
    for i=1:length(x)-1
        b = x(i+1)-x(i); h = airfoil_points(i+1)+airfoil_points(i); A = b*h;
        xA = (x(i+1)+x(i))/2;
        I = (h*b^3)/12;
        I_cm = I_cm + rho_r*(I + A*(xcm_r-xA)^2);
    end

    % Substraction of slot
    b = 0.03; h = 0.005; A_slot = b*h;
    xA = 0.05;
    I = (h*b^3)/12;    
    I_cm = I_cm - rho_r*(I + A_slot*(xcm_r-xA)^2);

    Icm_r = I_cm;
end
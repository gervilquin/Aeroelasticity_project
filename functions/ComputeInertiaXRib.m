function Ix_r = ComputeInertiaXRib(rho_r,t,c)
    I_cm = 0;

    % Computation without slot
    [x,airfoil_points] = CoordinatesSymmetricAirfoil(t,0.0075,0.095,200,c,'TE');
    for i=1:length(x)-1
        b = x(i+1)-x(i); h = airfoil_points(i+1)+airfoil_points(i);
        I = (b*h^3)/12;
        I_cm = I_cm + rho_r*I;
    end

    % Substraction of slot
    b = 0.03; h = 0.005;
    I = (b*h^3)/12;    
    I_cm = I_cm - rho_r*I;

    Ix_r = I_cm;
end
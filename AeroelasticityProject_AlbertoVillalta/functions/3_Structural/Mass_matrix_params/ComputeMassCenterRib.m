function xcm_r = ComputeMassCenterRib(A_r,t,c)
     x_cm = 0;
    
    % Computation without slot
    [x,airfoil_points] = CoordinatesSymmetricAirfoil(t,0.0075,0.095,200,c,'TE');
    for i=1:length(x)-1
        b = x(i+1)-x(i); h = airfoil_points(i+1)+airfoil_points(i); A = b*h;
        xA = (x(i+1)+x(i))/2;
        x_cm = x_cm + A*xA;
    end

    % Substraction of slot
    b = 0.03; h = 0.005; A_slot = b*h;
    xA = 0.05;
    x_cm = x_cm - A_slot*xA;

    xcm_r = x_cm/A_r;
end
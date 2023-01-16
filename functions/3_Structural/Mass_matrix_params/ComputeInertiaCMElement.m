function [I_cm] = ComputeInertiaCMElement(x_cm,material,t,c)
    I_cm = 0;
    
    % Leading edge
    [x,airfoil_points] = CoordinatesSymmetricAirfoil(t,0,0.0075,100,c,'LE');
    for i=1:length(x)-1
        b = x(i+1)-x(i); h = airfoil_points(i+1)+airfoil_points(i); A = b*h;
        xA = (x(i+1)+x(i))/2;
        I = (h*b^3)/12;
        I_cm = I_cm + material.Nylon.rho*(I + A*(x_cm-xA)^2);
    end

    % Plate
    b = 0.060; h = 0.00225; A = h*b;
    xA = 0.05;
    I = (h*b^3)/12;
    I_cm = I_cm + material.Al.rho*(I + A*(x_cm-xA)^2);

    % Trailing edge
    [x,airfoil_points] = CoordinatesSymmetricAirfoil(t,0.095,0.100,100,c,'TE');
    for i=1:length(x)-1
        b = x(i+1)-x(i); h = airfoil_points(i+1)+airfoil_points(i); A = b*h;
        xA = (x(i+1)+x(i))/2;
        I = (h*b^3)/12;
        I_cm = I_cm + material.Nylon.rho*(I + A*(x_cm-xA)^2);
    end

end
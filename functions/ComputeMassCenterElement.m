function [x_cm] = ComputeMassCenterElement(material,rhoA,t,c)
    x_cm = 0;

    % Leading edge
    [x,airfoil_points] = CoordinatesSymmetricAirfoil(t,0,0.0075,100,c,'LE');
    for i=1:length(x)-1
        b = x(i+1)-x(i); h = airfoil_points(i+1)+airfoil_points(i); A = b*h;
        xA = (x(i+1)+x(i))/2;
        x_cm = x_cm + material.Nylon.rho*A*xA;
    end

    % Plate
    b = 0.060; h = 0.00225; A = b*h;
    xA = 0.05;
    x_cm = x_cm + material.Al.rho*A*xA;

    % Trailing edge
    [x,airfoil_points] = CoordinatesSymmetricAirfoil(t,0.095,0.100,100,c,'TE');
    for i=1:length(x)-1
        b = x(i+1)-x(i); h = airfoil_points(i+1)+airfoil_points(i); A = b*h;
        xA = (x(i+1)+x(i))/2;
        x_cm = x_cm + material.Nylon.rho*A*xA;
    end
    
    x_cm = x_cm/rhoA;
end
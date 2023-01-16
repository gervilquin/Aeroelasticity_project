function [rhoA] = ComputeLinearDensityElement(material,t,c)
    rhoA = 0;

    % Leading edge
    [x,airfoil_points] = CoordinatesSymmetricAirfoil(t,0,0.0075,100,c,'LE');
    for i=1:length(x)-1
        b = x(i+1)-x(i); h = airfoil_points(i+1)+airfoil_points(i);
        A = b*h;
        rhoA = rhoA + material.Al.rho*A;
    end

    % Plate
    b = 0.060; h = 0.0025;
    A = b*h;
    rhoA = rhoA + material.Al.rho*A;

    % Trailing edge
    [x,airfoil_points] = CoordinatesSymmetricAirfoil(t,0.095,0.100,100,c,'TE');
    for i=1:length(x)-1
        b = x(i+1)-x(i); h = airfoil_points(i+1)+airfoil_points(i);
        A = b*h;
        rhoA = rhoA + material.Al.rho*A;
    end
end
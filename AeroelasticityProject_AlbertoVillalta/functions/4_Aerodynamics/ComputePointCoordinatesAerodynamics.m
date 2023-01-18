function p_coor = ComputePointCoordinatesAerodynamics(y_el,c,x_pos)
    % Initialisation
    nel = length(y_el)-1;
    p_coor = zeros(nel,3);

    % Coordinates computation
    for i=1:nel
        p_coor(i,1) = c*x_pos;               % x coordinate
        p_coor(i,2) = (y_el(i+1)+y_el(i))/2; % y coordinate
        p_coor(i,3) = 0;                     % z coordinate
        
    end
end
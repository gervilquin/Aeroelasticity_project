function p_coor = compute_aero_point(y_vector,chord_vector,x_pos)
    
    p_coor = zeros(length(y_vector)-1,3);


    for i=1:length(p_coor(:,1))
        p_coor(i,1) = chord_vector(i)*x_pos; % x coordinate
        p_coor(i,2) = (y_vector(i+1)-y_vector(i))/2 + y_vector(i); % y coordinate
        p_coor(i,3) = 0;
        
    end

end
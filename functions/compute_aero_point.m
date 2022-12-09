function [p_x,p_y] = compute_aero_point(y_vector,chord_vector,x_pos)
    
    p_y = zeros(1,length(y_vector)-1);
    p_x = zeros(1,length(y_vector)-1);

    for i=1:length(p_y)
        p_y(i) = (y_vector(i+1)-y_vector(i))/2 + y_vector(i);
        p_x(i) = chord_vector(i)*x_pos;
    end

end
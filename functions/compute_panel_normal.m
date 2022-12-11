function normal_vector = compute_panel_normal(alpha_aero)

    normal_vector = zeros(lenght(alpha_aero),2);
    for i=1:lenght(normal_vector(:,1))
        normal_vector(i,1) = cos(alpha_aero(i));
        normal_vector(i,2) = sin(alpha_aero(i)); 
    end

end
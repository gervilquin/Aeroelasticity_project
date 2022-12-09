function S = compute_element_surface(chord_vector,y_vector)
    S = zeros(1,length(y_vector)-1);
    for i=1:length(S)
        S(i) = (chord_vector(i+1)+chord_vector(i))/2 * (y_vector(i+1)-y_vector(i));
    end    
    
end
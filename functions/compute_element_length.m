function S = compute_element_length(y_vector)
    S = zeros(1,length(y_vector)-1);
    for i=1:length(S)
        S(i) = (y_vector(i+1)-y_vector(i));
    end   
    S = diag(S);%*1e-3;
    
end
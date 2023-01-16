function lambda_p = solve_system(D_prime)

    [V,D]=eigs(D_prime,eye(length(D_prime)),length(D_prime),'lr'); % lm?
    
    for k=1:size(V,2)
        lambda_p(k)=D(k,k);
    end
end
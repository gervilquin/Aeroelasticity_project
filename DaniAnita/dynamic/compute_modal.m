function [Phi,lambda] = compute_modal(K,M,modes,Ndof,If)

    K = (K + K')/2;
    M = (M + M')/2;
    
    [V,D] = eigs(K(If,If),M(If,If),modes,'sm');
    
    Phi = sparse(Ndof,modes);
    lambda = zeros(1,modes);
    
    for k = 1:size(V,2)
        Phi(If,k) = V(:,k)/sqrt((V(:,k)')*(M(If,If)*V(:,k)));
        lambda(1,k) = D(k,k);
    end
end
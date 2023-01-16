function [U] = reduced_order(K,M,modes,Ndof,If,f)
    [Phi,w2] = compute_modal(K,M,modes,Ndof,If);
    Im = 1:modes;
    
    U = zeros(Ndof,size(f,2));
    alpha = zeros(length(Im),size(f,2));
    
    for k = 1:size(f,2)
        for j = 1:length(Im)
            alpha(j,k) = (Phi(:,Im(j))')*f(:,k)/(w2(1,j));
            U(:,k) = U(:,k) + Phi(:,Im(j))*alpha(j,k);
        end
    end
end
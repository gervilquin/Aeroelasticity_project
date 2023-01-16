function L = compute_L(rho_inf,u_inf,S,alpha,A)
    S = diag(S);
    L = -rho_inf.*u_inf^2.*S*(A\alpha);
end
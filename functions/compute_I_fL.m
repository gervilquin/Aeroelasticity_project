function I_fL = compute_I_fL(ndof,nsec,e)

    %I_fL = sparse(ndof,nsec);
    I_fL = zeros(ndof,nsec);

    N_theta = [0.5 0.5];
    N_w = [0.5 0.5];
    N_gamma = [0.25 0.25];

    for i = 1:nsec 
        row = (i-1)*3;
        I_fL(row+1,i) = -e*N_theta(1);
        I_fL(row+2,i) = N_w(1);
        I_fL(row+3,i) = N_gamma(1);
        I_fL(row+4,i) = -e*N_theta(2);
        I_fL(row+5,i) = N_w(2);
        I_fL(row+6,i) = N_gamma(1);
    end
    
end


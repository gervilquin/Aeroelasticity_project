function I_fL = compute_I_fL(e)
    
    N_theta = [0.5 0.5];
    N_w = [0.5 0.5];
    N_gamma = [0.25 0.25];

    I_fL = [-e*N_theta(1) - N_theta(1);
            N_w(1) 0;
            N_gamma(1) 0;
            -e*N_theta(2) - N_theta(2);
            N_w(2) 0;
            N_gamma(2) 0];

end


function [I_alpha_u_0,I_alpha_u_1,I_alpha_u_2] = compute_I_au(U_inf,n,nsec,Tn)

   I_alpha_u_0 = sparse(nsec,3*n);
   I_alpha_u_1 = sparse(nsec,3*n);
   I_alpha_u_2 = sparse(nsec,3*n);

    for i = 1:nsec
        I = zeros(2*3,1); % (nnodes_el*nDOFs, 1)

        for k = 1:3 % 3 DOFs
            I(k,1) = 3*(Tn(i,1)-1)+k;
            I(k+3,1) = 3*(Tn(i,2)-1)+k;
        end

        I_alpha_u_0(i,I) = I_alpha_u_0(i,I) + [0.5 0 0 0.5 0 0];
        I_alpha_u_1(i,I) = I_alpha_u_1(i,I) -1/U_inf*[0 0.5 0 0 0.5 0];
        I_alpha_u_2(i,I) = I_alpha_u_2(i,I);
    end

end
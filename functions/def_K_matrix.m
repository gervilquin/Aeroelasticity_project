function K = def_K_matrix(y_nodal,Neset,Nsec,Nnod,Tn,EI,GJ)

    % Define a sparse K matrix
    %K = sparse(3*Nnod,3*Nnod);
    K = zeros(3*Nnod,3*Nnod);

    % Start the loop to fill the matrix
    for i = 1:Neset*Nsec

        I = zeros(2*3,1); % (nnodes_el*nDOFs, 1)
        l = y_nodal(Tn(i,2)) - y_nodal(Tn(i,1)); % element length

        for k = 1:3 % 3 DOFs
            I(k,1) = 3*(Tn(i,1)-1)+k;
            I(k+3,1) = 3*(Tn(i,2)-1)+k;
        end

        % Torsion stiffness matrix
        Kt = (GJ/l)*[ 1 0 0 -1 0 0;
                      0 0 0  0 0 0;
                      0 0 0  0 0 0;
                     -1 0 0  1 0 0;
                      0 0 0  0 0 0;
                      0 0 0  0 0 0];

        % Bending stifness matrix
        Kb = (EI/l^3)*[0  0    0   0   0   0   ;
                       0  12  6*l  0  -12  6*l ; 
                       0 6*l 4*l^2 0 -6*l 2*l^2;
                       0  0    0   0   0   0   ;
                       0 -12 -6*l  0   12  -6*l;
                       0 6*l 2*l^2 0 -6*l 4*l^2];
        
        % Add the element matrix to the global matrix 
        K(I,I) = K(I,I) + Kt + Kb;
    end

end
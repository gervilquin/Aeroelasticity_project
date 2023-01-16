function K = ComputeKmatrix(y_el,Tn,EI,GJ)
% Definition and computation of structural model stiffness matrix

    % Initialisation
    nodeDOFs = 3;
    nel = size(Tn,1);
    nnod = length(y_el);
    K = sparse(nodeDOFs*nnod,nodeDOFs*nnod);
 
    % Matrix computation
    for i = 1:nel
        l = y_el(Tn(i,2)) - y_el(Tn(i,1)); % Element's length

        % DOF determination for allocation
        I = zeros(size(Tn,2)*nodeDOFs,1); % [nnodes_el*nDOFs, 1]
        for k = 1:nodeDOFs
            I(k,1) = 3*(Tn(i,1)-1)+k;
            I(k+3,1) = 3*(Tn(i,2)-1)+k;
        end

        % Torsional rigidity matrix
        Kt = (GJ/l)*[ 1 0 0 -1 0 0;
                      0 0 0  0 0 0;
                      0 0 0  0 0 0;
                     -1 0 0  1 0 0;
                      0 0 0  0 0 0;
                      0 0 0  0 0 0];

        % Flexural rigidity matrix
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
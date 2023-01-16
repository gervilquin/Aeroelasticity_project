function I_fL = ComputeForcesCoupling(y_el,Tn,e)
    % Initialisation
    nodeDOFs = 3;
    nnodes = length(y_el);
    nDOFs = nodeDOFs*nnodes;
    nel = nnodes - 1;
    I_fL = sparse(nDOFs,nel);

    % Shape function values at element's center
    N_theta = [0.5 0.5];
    N_w = [0.5 0.5];
    N_gamma = [0.25 -0.25];

    % Computation
    for i = 1:nel 
        % DOFs of the element
        l = y_el(Tn(i,2)) - y_el(Tn(i,1)); % Element's length
        I = zeros(size(Tn,2)*nodeDOFs,1); % [nnodes_el*nDOFs, 1]
        for k = 1:nodeDOFs
            I(k,1) = 3*(Tn(i,1)-1)+k;
            I(k+3,1) = 3*(Tn(i,2)-1)+k;
        end
        % Coupling values
        I_fL(I,i) = I_fL(I,i) + [-e*N_theta(1); N_w(1); N_gamma(1)*(l/2); -e*N_theta(2); N_w(2); N_gamma(2)*(l/2)];
    end
    
end


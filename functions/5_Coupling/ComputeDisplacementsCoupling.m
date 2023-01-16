function [I_alpha_u_1, I_alpha_u_2, I_alpha_u_3] = ComputeDisplacementsCoupling(y_el,Tn,U_inf)
    % Initialisation
    nodeDOFs = 3;
    nnodes = length(y_el);
    nDOFs = nodeDOFs*nnodes;
    nel = nnodes - 1;
    I_alpha_u_1 = sparse(nel,nDOFs);
    I_alpha_u_2 = sparse(nel,nDOFs);
    I_alpha_u_3 = sparse(nel,nDOFs);

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
        % Displacements
        I_alpha_u_1(i,I) = I_alpha_u_1(i,I) + [N_theta(1), 0, 0, N_theta(2), 0, 0];
        % Velocities
        I_alpha_u_2(i,I) = I_alpha_u_2(i,I) + (-1/U_inf) * [0, N_w(1), N_gamma(1)*(l/2), 0, N_w(2), N_gamma(2)*(l/2)];
        % Accelerations
        I_alpha_u_3(i,I) = I_alpha_u_3(i,I) + zeros(1,6);

    end
end
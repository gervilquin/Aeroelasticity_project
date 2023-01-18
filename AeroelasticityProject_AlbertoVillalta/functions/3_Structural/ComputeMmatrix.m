function M = ComputeMmatrix(y_el,Tn,Ts,material,t,c,x_sc,RibExist)
% Definition and computation of structural model mass matrix

    % Initialisation
    nodeDOFs = 3;
    nel = size(Tn,1);
    nnod = length(y_el);
    M = zeros(nodeDOFs*nnod,nodeDOFs*nnod);

    % Matrix computation (Beam elements)
    [rhoA_e] = ComputeLinearDensityElement(material,t,c);    % Linear density (kg/m)
    [xcm_e] = ComputeMassCenterElement(material,rhoA_e,t,c); % Location of center of mass (m)
    [Icm_e] = ComputeInertiaCMElement(xcm_e,material,t,c);   % Section's moment of inertia at center of mass (kg·m)

    for i = 1:nel
        l = y_el(Tn(i,2)) - y_el(Tn(i,1)); % Element's length

        % DOF determination for allocation
        I = zeros(size(Tn,2)*nodeDOFs,1); % [nnodes_el*nDOFs, 1]
        for k = 1:nodeDOFs
            I(k,1) = nodeDOFs*(Tn(i,1)-1)+k;
            I(k+3,1) = nodeDOFs*(Tn(i,2)-1)+k;
        end        
    
        % Torsion matrix
        Mt = (Icm_e*l/6)*[2 0 0 1 0 0;
                          0 0 0 0 0 0;
                          0 0 0 0 0 0;
                          1 0 0 2 0 0;
                          0 0 0 0 0 0;
                          0 0 0 0 0 0];
        % Bending matrix
        Mb = (rhoA_e*l/420)*[0,   0,       0,   0,    0,       0 ;
                             0,  156,    22*l,  0,   54,    -13*l;
                             0,  22*l,  4*l^2,  0,  13*l,  -3*l^2;
                             0,   0,       0,   0,    0,       0 ;
                             0,   54,    13*l,  0,   156,   -22*l;
                             0, -13*l, -3*l^2,  0, -22*l,   4*l^2];

        d = xcm_e - x_sc;%*c;
        d_m = [ 1  0  0  0  0  0;
               -d  1  0  0  0  0;
                0  0  1  0  0  0;
                0  0  0  1  0  0;
                0  0  0 -d  1  0;
                0  0  0  0  0  1];

        M(I,I) = M(I,I) + d_m'*(Mt + Mb)*d_m;
    end

    % Matrix computation (Punctual masses)
    sec_ribs = find(RibExist == 1);  % Sections with ribs
    Nribs = length(sec_ribs);        % Number of ribs
    h_r = 0.004;                     % Thickness of the ribs (m)
    rho_r = material.Nylon.rho;      % Rib density (made of nylon)

    [A_r]   = ComputeAreaRib(t,c);                   % Area of the rib (m^2)
    [xcm_r] = ComputeMassCenterRib(A_r,t,c);         % Rib's center of mass (m)
    [Icm_r] = ComputeInertiaCMRib(rho_r,xcm_r,t,c);  % Rib's moment of inertia at center of mass (kg·m)
    [Ix_r]  = ComputeInertiaXRib(rho_r,t,c);         % Rib's moment of inertia in x-axis (kg·m)

    Mr = h_r*[Icm_r+rho_r*A_r*(xcm_r-x_sc)^2, -rho_r*A_r*(xcm_r-x_sc),   0;
               -rho_r*A_r*(xcm_r-x_sc),              rho_r*A_r,         0;
                          0,                           0,           Ix_r];
    for j = 1:Nribs
        % DOF determination for allocation
        I = zeros(1,nodeDOFs);
        for k =1:nodeDOFs
            I(k) = nodeDOFs*(Tn(Ts(sec_ribs(j),1),2) -1) + k;
        end
        M(I,I) = M(I,I) + Mr;
    end
    

end
function M = def_M_matrix(y_nodal,Neset,Nsec,Nnod,Tn,Ts,material,naca_2dt,chord,x_sc)

    % Initialize the mass matrix
    %M = sparse(3*Nnod,3*Nnod);
    M = zeros(3*Nnod,3*Nnod);

    % Start loop for LE, TE and Aluminium plate
    for i = 1:Neset*Nsec

        I = zeros(2*3,1); % (nnodes_el*nDOFs, 1)
        l = y_nodal(Tn(i,2)) - y_nodal(Tn(i,1)); % element length

        for k = 1:3 % 3 DOFs
            I(k,1) = 3*(Tn(i,1)-1)+k;
            I(k+3,1) = 3*(Tn(i,2)-1)+k;
        end

        [rhoA] = compute_avg_density(material,naca_2dt,chord);
        [x_cm] = compute_mass_center(material,rhoA,naca_2dt,chord);
        [I_cm] = compute_inertia_mass_center(x_cm,material,naca_2dt,chord);
    
        % Torsion matrix
        Mt = (I_cm*l/6)*[2 0 0 1 0 0;
                        0 0 0 0 0 0;
                        0 0 0 0 0 0;
                        1 0 0 2 0 0;
                        0 0 0 0 0 0;
                        0 0 0 0 0 0];
        % Bending matrix
        Mb = (rhoA*l/420)*[0   0       0   0    0       0 ;
                           0  156    22*l  0   54    -13*l;
                           0  22*l  4*l^2  0  13*l  -3*l^2;
                           0   0       0   0    0       0 ;
                           0   54    13*l  0   156   -22*l;
                           0 -13*l -3*l^2  0 -22*l   4*l^2];

        d = x_cm - x_sc;
        d_m = [ 1  0  0  0  0  0;
             -d  1  0  0  0  0;
              0  0  1  0  0  0;
              0  0  0  1  0  0;
              0  0  0 -d  1  0;
              0  0  0  0  0  1];
        
        M(I,I) = M(I,I) + d_m'*(Mt + Mb)*d;

    end

    % Ribs
    Nribs = length(Ts(:,1));
    h_r = 0.004; %mm thickness of the ribs
    rho_r = material.Nylon.rho;

    % Start loop for ribs
    for j = 1:Nribs

        I = zeros(1,3);
        for k =1:3
            I(k) = 3*(Tn(Ts(j,2),2) -1) + k;
        end

        A_r = compute_rib_area(naca_2dt,chord);
        xcm_r= compute_rib_xcm(A_r,naca_2dt,chord); 
        Icm_r= compute_rib_Icm(rho_r,xcm_r,naca_2dt,chord);
        Ix_r = compute_rib_Ix(rho_r,naca_2dt,chord);

        Mr = h_r*[Icm_r+rho_r*A_r*(xcm_r-x_sc)^2,     -rho_r*A_r*(xcm_r-x_sc),    0;
                  -rho_r*A_r*(xcm_r-x_sc)           rho_r*A_r                   0;
                    0                               0                           Ix_r];

        M(I,I) = M(I,I) + Mr;
    end
    
disp("Slender Beam:")
disp(["rhoA",num2str(rhoA)])
disp(["x_cm",num2str(x_cm)])
disp(["I_cm",num2str(I_cm)])
 
disp("Rib")
disp(["rhoA",num2str(rho_r*A_r)])
disp(["xcm_r",num2str(xcm_r)])
disp(["Icm_r",num2str(Icm_r)])
disp(["Ix_r",num2str(Ix_r)])


end
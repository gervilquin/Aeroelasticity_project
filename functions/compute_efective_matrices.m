function [Meff,Ceff,Keff] = compute_efective_matrices(M,K,U_inf,Parameters)

    
    S = Parameters.S;
    rho_inf = Parameters.rho;
    Nnod = Parameters.Nnod;
    Nel = Parameters.Nel;
    Tn = Parameters.Tn ;
    If = Parameters.If;
    I_fL = Parameters.I_fl;
    A_aero = Parameters.A_aero;
    
    % Compute Aerodynamic matrices
    S_aero = -U_inf^2*rho_inf*S;

    % Compute coupling matrices
    [I_au_0,I_au_1,I_au_2] = compute_I_au(U_inf,Nnod,Nel,Tn);
    %I_fL = compute_I_fL(Nnod,Nel,x_sc-x_ac);

%     I_au_0 = I_au_0(:,If);
%     I_au_1 = I_au_1(:,If);
%     I_au_2 = I_au_2(:,If);

    % Compute aerodynamic matrices
    M_a = I_fL*(S_aero*inv(A_aero))*I_au_2;
    C_a = I_fL*(S_aero*inv(A_aero))*I_au_1;
    K_a = I_fL*(S_aero*inv(A_aero))*I_au_0;

    % Compute efective matrices
    Meff = M + M_a; 
    Ceff = C_a;
    Keff = K + K_a;
    



end
function [Meff,Ceff,Keff] = ComputeEffectiveMatrix(K,M,y_el,Tn,U_inf,rho_inf,A_aero,S,I_fL)
    % Compute velocity dependant matrices
    % Compute aerodynamic matrices
    S_aero = (-U_inf^2)*rho_inf*S;
    % Compute coupling matrices
    [I_au_0, I_au_1, I_au_2] = ComputeDisplacementsCoupling(y_el,Tn,U_inf);  

    % Compute aerodynamic matrices
    M_a = -I_fL*(S_aero*inv(A_aero))*I_au_2;
    C_a = -I_fL*(S_aero*inv(A_aero))*I_au_1;
    K_a = -I_fL*(S_aero*inv(A_aero))*I_au_0;

    % Compute effective matrices
    Meff = M + M_a; 
    Ceff = C_a;
    Keff = K + K_a;
end
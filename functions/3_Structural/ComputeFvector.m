function [F,L] = ComputeFvector(nel,U_inf,rho_inf,AoA_deg,S,A_aero,I_fL)
    % Compute force vector
    AoA = deg2rad(AoA_deg);
    alpha = AoA*ones(nel,1);
    S_aero = -(U_inf^2)*rho_inf*S;
    L = S_aero*inv(A_aero)*alpha;
    F= I_fL*L;
end
%% Compute S aerodynamic matrix
%   input
% u_inf    -> freestream velocity
% rho_inf  -> freestream density
% Sf       -> Surface of each partition   
%   output
%         by comparison of slide 22 project with slide 1 S06 I
%         understand that S has this value
%         L = rhoinf*Uinf*S[i]*Gamma[i] <-> {L} = [S(Uinf)]{Gamma}
% S        -> S matrix

% Authors: Ana María Díaz, Daniel Regener
% Advanced Aeroelasticity, ESEIAAT, 31 dec 2022
function S = compute_S(rho_inf,l)
    S_aux = zeros(length(l),length(l));
    for i = 1:length(l)
        S_aux(i,i) = rho_inf*l(i);
    end
    S = @(u_inf) S_aux.*u_inf;
end
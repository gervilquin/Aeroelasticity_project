%% Compute D' matrix for the stability analysis
%   input
% C      -> reduced system effective stifness matrix
% C      -> reduced system effective damping matrix
% M      -> reduced system effective mass matrix
%   output
%         from S06 - slide 13, we can calulate D'
% D_prime -> D' matrix

% Authors: Ana María Díaz, Daniel Regener
% Advanced Aeroelasticity, ESEIAAT, 31 dec 2022
function D_prime = compute_D_prime(K,C,M)
    D_prime = [K\C       ,K\M;...
               -eye(size(K,1)),zeros(size(K,1))];
end
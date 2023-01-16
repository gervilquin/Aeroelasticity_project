%% Mass normalization of eigenvectors of a system
%   input
% K     -> Stifness matrix
% M     -> Mass matrix
%   output
% phi   -> Mass normalized eigenvectors
%   Notes
% As proposed in S06 - slide 12
% Authors: Ana María Díaz, Daniel Regener
% Advanced Aeroelasticity, ESEIAAT, 28 dec 2022

function phi = mass_normalization(phi,M)
    phi = sqrt(phi'*M*phi)^(-1)*phi;
end
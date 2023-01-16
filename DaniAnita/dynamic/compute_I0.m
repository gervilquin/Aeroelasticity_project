%% Create I0 matrix for apying BC
%   input
% prescribed_dofs  -> int, number of prescribed dofs
% ndofs            -> int, total number of dofs
%   output
% I0                -> As proposed in S06 - slide 2ç
%   example of call
% I0 = compute_I0(3,100)
% Authors: Ana María Díaz, Daniel Regener
% Advanced Aeroelasticity, ESEIAAT, 31 dec 2022
function I0 = compute_I0(prescribed_dofs,ndofs)
    I0 = zeros(ndofs,ndofs-prescribed_dofs);
    I0(prescribed_dofs+1:end,1:end) = eye(ndofs-prescribed_dofs);
end
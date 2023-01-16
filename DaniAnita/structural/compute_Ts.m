%% Create nodal connectivity matrix
%   input
% ns           -> number of nodes in each subset (int vector)
%   output
% Ts           -> nodal connectivity matrix

% Authors: Ana María Díaz, Daniel Regener
% Advanced Aeroelasticity, ESEIAAT, 28 dec 2022
function Ts = compute_Ts(ns)
    count = 1;
    for i = 1:size(ns,2)
        for j = 1:ns(i)-1
           Ts(i,j) = count;
           count = count + 1;
        end
    end
end
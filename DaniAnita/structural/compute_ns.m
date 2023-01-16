%% Compute ns vector
%   input
% nodes_subset -> nodes per subset (int)
% nodes_rib    -> nodes in ribs    (int)
% dead_ribs    -> ribs that has no node_rib after them (int)
%   output
% ns           -> number of nodes in each subset (int vector)

% Authors: Ana María Díaz, Daniel Regener
% Advanced Aeroelasticity, ESEIAAT, 28 dec 2022
function ns = compute_ns(nodes_subset,nodes_rib,dead_ribs)
    for e=1:nodes_rib-1
        if e >= (nodes_rib-1)-(dead_ribs-1)
            ns(e) = 2;
        else
            ns(e) = nodes_subset; 
        end
    end
end
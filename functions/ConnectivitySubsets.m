function [Ts] = ConnectivitySubsets(Tn, nset, neset)
% Connectivity table of the subsets (2 nodes per subset)
% Ts [#subsets, #nodes per subset considered]
    Ts = zeros(nset,2);
    for i=1:nset
        Ts(i,:) = [Tn((i-1)*neset+1,1), Tn(i*neset,2)];
    end
end
function [Ts] = ConnectivitySubsets(Tn, nset, neset)
% Connectivity table of the subsets (2 nodes per subset)
% Ts [#subsets, #nodes per subset considered]
    Ts = zeros(nset,2);
    for i=1:nset-1
        Ts(i,:) = [Tn(neset*(i-1)+1,1), Tn(i*neset + 1,1)];
    end
    Ts(end,:) = [Tn(neset*(nset-1)+1,:)];
end
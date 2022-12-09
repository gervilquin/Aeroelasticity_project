function [Tn] =ConnectivityElements(nel)
% Connectivity table of the elements (2 nodes per element)
% Tn [#elements, #nodes per element]
    Tn = zeros(nel,2);
    for i=1:nel
        Tn(i,:) = [i, i+1];
    end
end
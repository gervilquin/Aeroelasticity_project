function S = ComputeElementSize(y_el)
    % Initialisation
    nel = length(y_el)-1;
    S = zeros(1,nel);

    % Computation
    for i=1:nel
        S(i) = (y_el(i+1)-y_el(i));
    end   
    S = diag(S);
end
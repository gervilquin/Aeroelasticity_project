%% Create element connectivity matrix
%   input
% ns           -> number of nodes in each subset (int vector)
%   output
% Tn           -> element connectivity matrix

% Authors: Ana María Díaz, Daniel Regener
% Advanced Aeroelasticity, ESEIAAT, 28 dec 2022

function [aux, Tn] = compute_Tn(ns)
    count = 1;
    for e=1:size(ns,2)
       for i = 1:ns(e)
           aux(e,i)=  count;
           if i == ns(e), break; end
           count = count +1;
            
       end
    end
    
    count = 1;
    for i=1:max(max(aux))-1
       for j = 1:2
           Tn(i,j)=  count;
           if j == 2, break; end
           count = count +1;
            
       end
    end
end

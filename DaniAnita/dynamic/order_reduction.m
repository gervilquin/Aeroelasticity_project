%% Compute eigenvectors and their frecuecies
%   input
% K     -> Stifness matrix
% M     -> Mass matrix
% Modes -> Number of modes desired for phi and w2
%   output
% phi   -> First eigenvectors to Modes
% w2    -> First eigenfrequencies to Modes
%   Notes
% Notes: This function assumes that boundary conditions are already applied
% Authors: Ana María Díaz, Daniel Regener
% Advanced Aeroelasticity, ESEIAAT, 31 dec 2022
function [phi,w2]=order_reduction(K,M,Modes)
    
    K=(K+K')/2;
    M=(M+M')/2;

    [V,D]=eigs(K,M,Modes,'sm');
    
    phi=zeros(size(K,1),Modes);
    w2=zeros(1,Modes);
    
    for k=1:size(V,2)
        phi(:,k)=V(:,k)/sqrt(V(:,k)'*M*V(:,k));
        w2(k)=D(k,k);
    end

end
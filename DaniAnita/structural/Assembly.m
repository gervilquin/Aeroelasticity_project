%% Assemble structural matrices
%   input
% y           -> position of the nodes
% l           -> length of each element
% GJ          -> tornsional stifness
% EI          -> bending stifness
% rho_al      -> density of the core
% rho_al      -> density of skins and front and rear spar
% Tn          -> element connectivity matrix
% ns          -> number of nodes in i subset
% xsc         -> shear center position
% Ts          -> nodal connectivity matrix
% myNACA      -> fucntion of x/c for profile of the wing, MUST BE SIMMETRIC
% hr          -> thickness of the ribs
%   output
% DOFs        -> total number of degrees of freedom
% K           -> global stifness matrix
% M           -> global mass matrix
% Authors: Ana María Díaz, Daniel Regener
% Advanced Aeroelasticity, ESEIAAT, 28 dec 2022
function [DOFs, K, M]  = Assembly(y, l,GJ,EI,rho_al,rho_po,Tn,ns,xsc,Ts,myNACA,hr)
    
    % calculate the total number of degrees of freedom
    DOFs = 3*(length(y)); % total numberze of DOFs

    % calculate inertia of section with core, front spar and rear spar
    [Icm, xcm,rhoA] = calcInertia(rho_al,rho_po,myNACA);
    
    % initialize stifness and mass matrices
    K = sparse(DOFs,DOFs);
    M = sparse(DOFs,DOFs);
    for i=1:length(l)
        for k=1:3
            % conversion from local DOF to global DOF
            I(k)   = 3*(Tn(i,1)-1) + k; 
            I(k+3) = 3*(Tn(i,2)-1) + k; 
        end
        % compute the element matrix for the i element
        [Ki,Mi] = elementMatrices(GJ,EI,l(i),Icm,xcm,xsc,rhoA);
        % assemble them into the global matrice
        K(I,I) = K(I,I) + Ki;
        M(I,I) = M(I,I) + Mi;
    end
    
    check_simmetry(M,'M without ribs');
    % Ribs contribution
    Mr = calcMr(rho_po,xsc,hr);
    nodesWithRib=ns(1);
    for e = 2:size(ns,2)-1
       nodesWithRib = [nodesWithRib, nodesWithRib(end)+(ns(e))]; 
    end
    
    for j = 1:size(Ts,1)-1 % the last node has no rib!!
        for k=1:3
            Ir(k) = 3*(Tn(Ts(j,1),2)-1)+k;
        end
        M(Ir,Ir) = M(Ir,Ir) + Mr;
    end
    

end
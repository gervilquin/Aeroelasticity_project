%% Assemble structural matrices
%   input
% l           -> length of each element
% GJ          -> tornsional stifness
% EI          -> bending stifness
% Icm         -> inertia in the center of mass
% xc,         -> x position of the center of mass
% rhoA        -> density per surface of each section
%   output
% K           -> element stifness matrix
% M           -> element mass matrix
% Authors: Ana María Díaz, Daniel Regener
% Advanced Aeroelasticity, ESEIAAT, 28 dec 2022
function [K,M] = elementMatrices(GJ,EI,l,Icm,xcm,xsc,rhoA)
    % Element torsional stifness matrix, slide 16
    Kt = zeros(6);
    Kt(1,1) = 1;
    Kt(1,4)=-1;
    Kt(4,1)=-1;
    Kt(4,4)=1;
    Kt = (GJ/l).*Kt;
    
    % Element torsional bending matrix, slide 16
    Kb= [0, 0,          0,      0,      0,      0;...
        0,  12,         6*l,    0,      -12,    6*l;...
        0,  6*l,        4*l^2,  0,      -6*l,   2*l^2;...
        0,  0,          0,      0,      0,      0;...
        0,  -12,        -6*l,   0,      12,    -6*l;...
        0,  6*l,        2*l^2,  0,      -6*l,   4*l^2];
    Kb = (EI/l^3).*Kb;
    
    K = Kt + Kb;
    
    Mtp = zeros(6);
    
    Mtp(1,1)=2;
    Mtp(1,4)=1;
    Mtp(4,1)=1;
    Mtp(4,4)=2;
    Mtp = ((Icm*l)/6).*Mtp;
    if (Mtp-Mtp')~=0, error('Typo! Non symmetric Mtp'); end
    
    Mbp =[0, 0,          0,      0,      0,      0;...
        0,  156,        22*l,    0,      -54,    -13*l;...
        0,  22*l,        4*l^2,  0,      13*l,   -3*l^2;...
        0,  0,          0,      0,      0,      0;...
        0,  54,        13*l,   0,      156,    -22*l;...
        0,  -13*l,        -3*l^2,  0,      -22*l,   4*l^2];
    Mbp = Mbp.*((rhoA*l)/420);
    if (Mbp-Mbp')~=0, error('Typo! Non symmetric Mbp'); end
    
    d=eye(6);
    d(2,1)=-(xcm-xsc);
    d(5,4)=-(xcm-xsc);
    
    M = d'*(Mtp+Mbp)*d;
end
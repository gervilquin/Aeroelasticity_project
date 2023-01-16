%% calculate inertia of section with core, front spar and rear spar
%   input
% rho_al      -> density of the core
% rho_al      -> density of skins and front and rear spar
%   output
% Icm         -> inertia in the center of mass
% xc,         -> x position of the center of mass
% rhoA        -> density per surface of each section
% Authors: Ana María Díaz, Daniel Regener
% Advanced Aeroelasticity, ESEIAAT, 28 dec 2022
% A Paola le da 7.15e-4 ??
function [Icm, xcm,rhoA] = calcInertia(rho_al,rho_po,myNACA)
   % compute the surfaces of each subsection (front spar, core, rear spar)
   S(1) = 2.*integral(myNACA,0,7.5e-3);
   S(2) = 2.25e-3.*60e-3;
   S(3) = 2.*integral(myNACA,95e-3,100e-3);
   % add them and obtain rhoA
   rhoA = S(1)*rho_po+S(2)*rho_al+S(3)*rho_po;
   if rhoA<0, error("Negative rhoA = %f",rhoA); end
   % compute the center of mass (slide 17)
   myIntegral = @(x) 2.*x.*myNACA(x); %dx
   xcm = rho_po.*integral(myIntegral,0,7.5e-3);
   myIntegral = @(x) 2.*x.*1.125e-3; %dx
   xcm = xcm+rho_al.*integral(myIntegral,20e-3,80e-3);
   myIntegral = @(x) 2.*x.*myNACA(x);
   xcm = xcm+rho_po.*integral(myIntegral,95e-3,100e-3);
   xcm = (1/rhoA).*xcm;
   % compute inertia in the center of mass (slide 17)
   myIntegral = @(x) 2.*(x-xcm).^2.*myNACA(x);
   Icm = rho_po.*integral(myIntegral,0,7.5e-3);
   myIntegral = @(x) 2.*(x-xcm).^2.*1.125e-3;
   Icm = Icm+rho_al.*integral(myIntegral,20e-3,80e-3);
   myIntegral = @(x) 2.*(x-xcm).^2.*myNACA(x);
   Icm = Icm+rho_po.*integral(myIntegral,95e-3,100e-3);
   if (xcm<0)||(xcm>0.1), error("Out of section xcm = %f",xcm); end
end
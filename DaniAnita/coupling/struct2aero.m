function [I_alpha_u0, I_alpha_u1, I_alpha_u2] = struct2aero(y,aerodynamic_point,Ndofs,xsc,xac)
    nt1 = @(eta) 0.5*(1-eta);
    nt2 = @(eta) 0.5*(1+eta);
    nw1 = @(eta) 0.25*(2-3*eta+eta^3);
    dnw1= @(eta) 0.75*eta^2-0.75;
    nw2 = @(eta) 0.25*(2+3*eta-eta^3);
    dnw2= @(eta) 0.75-0.75*eta^2;
    ng1 = @(eta) 0.25*(1-eta-eta^2+eta^3);
    dng1= @(eta) 0.75*eta^2-0.5*eta+0.25;
    ng2 = @(eta) 0.25*(-1-eta+eta^2+eta^3);
    dng2= @(eta) 0.75*eta^2-0.5*eta-0.25;
    eta_f = @(yc,y,yp1) (2*yc-(y+yp1))/(yp1-y);

    I_star_u = zeros(Ndofs-3,Ndofs);
    for i=1:length(aerodynamic_point)
       eta = eta_f(aerodynamic_point(i,2),y(i),y(i+1));
       I_star_u(1+3*(i-1),1+3*(i-1)) = nt1(eta);
       I_star_u(1+3*(i-1),4+3*(i-1)) = nt2(eta);
       I_star_u(2+3*(i-1),2+3*(i-1)) = nw1(eta);
       I_star_u(2+3*(i-1),3+3*(i-1)) = ng1(eta);
       I_star_u(2+3*(i-1),5+3*(i-1)) = nw2(eta);
       I_star_u(2+3*(i-1),6+3*(i-1)) = ng2(eta);
       I_star_u(3+3*(i-1),2+3*(i-1)) = dnw1(eta);
       I_star_u(3+3*(i-1),3+3*(i-1)) = dng1(eta);
       I_star_u(3+3*(i-1),5+3*(i-1)) = dnw2(eta);
       I_star_u(3+3*(i-1),6+3*(i-1)) = dng2(eta);
    end
    
    f = 1/(xsc-xac);
    B = zeros(length(aerodynamic_point),Ndofs-3);
    count = 1;
     for i=1:length(aerodynamic_point)
       B(i,count)   = 1;
%        B(i,count+1) = -f;
       B(i,count+1) = 0;
       B(i,count+2) = 0;
       count = count+3;
     end
    I_alpha_u0_aux = B*I_star_u;
    I_alpha_u0 = @(U_inf) I_alpha_u0_aux;

    
    I_alpha_u1_aux = zeros(length(aerodynamic_point),Ndofs);
    for i=1:length(aerodynamic_point)
       eta = eta_f(aerodynamic_point(i,2),y(i),y(i+1));
       I_alpha_u1_aux(i,2+3*(i-1)) = nw1(eta);
       I_alpha_u1_aux(i,3+3*(i-1)) = ng1(eta);
       I_alpha_u1_aux(i,5+3*(i-1)) = nw2(eta);
       I_alpha_u1_aux(i,6+3*(i-1)) = ng2(eta);
    end
    I_alpha_u1 =  @(U_inf) -(1/U_inf).*I_alpha_u1_aux;

    I_alpha_u2 =  @(U_inf) zeros(size(I_alpha_u1_aux ));
end
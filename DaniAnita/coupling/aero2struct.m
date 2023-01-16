function [I_f_L] = aero2struct(y,aerodynamic_point,DOFs,xsc,xac)
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
    f     = - xac + xsc;
    I_f_L = zeros(DOFs,length(aerodynamic_point));
    
    
    for i = 1:length(aerodynamic_point)
        if i==1
            eta = 1; % creo k no importa, porque es el nodo restringido
%             I_f_L(1,1) = nt2(eta);
%             I_f_L(2,1) = nw2(eta);
              I_f_L(1,1) = 0;
              I_f_L(2,1) = 0;
              I_f_L(3,1) = 0;
%         %elseif i ==length(y)
%             % ????
        else
%             eta = eta_f(y(i),aerodynamic_point(i-1,2),aerodynamic_point(i,2));
%             I_f_L(1+3*(i-1),i-1) = nt1(eta)*f;
%             I_f_L(1+3*(i-1),i) = nt2(eta)*f;
%             I_f_L(2+3*(i-1),i-1) = nw1(eta);
%             I_f_L(2+3*(i-1),i) = nw2(eta);
            Delta_y_L(1) = abs(aerodynamic_point(i-1,2)-y(i));
            Delta_y_L(2) = abs(aerodynamic_point(i,2)-y(i));
            I_f_L(1+3*(i-1),i-1) = 0.5*f;
            I_f_L(1+3*(i-1),i) = 0.5*f;
            I_f_L(2+3*(i-1),i-1) = 0.5;
            I_f_L(2+3*(i-1),i) =0.5;
            I_f_L(3+3*(i-1),i-1) =  -0.5*Delta_y_L(1);
            I_f_L(3+3*(i-1),i)   =   0.5*Delta_y_L(2);

        end
        
    end

end
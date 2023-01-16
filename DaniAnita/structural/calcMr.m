function Mr = calcMr(rho,xsc,hr)
    h = 5e-3;
    t = 0.018; % Aqui hi havia un 0.18
    myNACA = @(x) 5.*t.*(0.2969.*sqrt(x/0.1)-0.1260.*(x/0.1)-0.3516.*(x/0.1).^2+...
       0.2843.*(x/0.1).^3-0.1015.*(x/0.1).^4); 
   
    % A = 2*integral(myNACA,0,100e-3);

    A = 2*integral(myNACA,7.5e-3,95e-3)-h*(65e-3-35e-3);

    myIntegral1 = @(x) 2.*x.*myNACA(x);
    myIntegral2 = @(x) x.*h;
    integrated = integral(myIntegral1,7.5e-3,95e-3)-integral(myIntegral2,35e-3,65e-3);
    xcmr = 1/(A).*integrated;
    
    myIntegral1 = @(x) 2.*rho.*(x-xcmr).^2.*myNACA(x); 
    myIntegral2 = @(x) rho.*(x-xcmr).^2.*h; 
    integrated = integral(myIntegral1,7.5e-3,95e-3)-integral(myIntegral2,35e-3,65e-3);
    Icmr = integrated;
    
    % He canviat la definició de Ixr pq es integral de: ((1/3)*y^3dx)
    myIntegral = @(x) (1/3)*myNACA(x).^3;
    Ixr = rho*(2*integral(myIntegral,7.5e-3,95e-3)-(1/3)*(95e-3-7.5e-3)*h^3);
    
    Mr(1,1) = Icmr+rho.*A.*(xcmr-xsc).^2;
    Mr(1,2) = -rho.*A.*(xcmr-xsc);
    Mr(1,3) = 0;
    Mr(2,1) = -rho.*A.*(xcmr-xsc);
    Mr(2,2) = rho.*A;
    Mr(2,3) = 0;
    Mr(3,1) = 0;
    Mr(3,2) = 0;
    Mr(3,3) = Ixr;
    
    Mr = hr.*Mr;
    
end


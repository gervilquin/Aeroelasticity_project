function [y,x] = compute_symmetric_airfoil(naca_2dt,x0,xf,n,length_airfoil)
    x = linspace(x0,xf,n)/length_airfoil;
    y = zeros(1,n);
    t = naca_2dt/100;
    for i=1:length(x)
        y(i) = 5*t*(0.2969*sqrt(x(i)) - 0.1260*x(i) - 0.3516*x(i)^2 + ...
                    + 0.2843*x(i)^3 - 0.1036*x(i)^4);
    end
    y = y*length_airfoil;
    x = x*length_airfoil;
end
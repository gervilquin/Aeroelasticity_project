function [x,y] = CoordinatesSymmetricAirfoil(t,x0,xf,n,c,mode)

    % To improve precision at leading edge high curvature
    if strcmp(mode,'LE')
        x = xf*(1-cos(linspace(0,pi/2,n)));
    elseif strcmp(mode,'TE')
        x = linspace(x0,xf,n)/c;
    end

    y = zeros(1,n);
    t = t/100; % Convert from percentage to decimal
    for i=1:length(x)
        y(i) = 5*t*(0.2969*sqrt(x(i)) - 0.1260*x(i) - 0.3516*x(i)^2 + ...
                    + 0.2843*x(i)^3 - 0.1036*x(i)^4);
    end
    y = y*c;
    x = x*c;

end
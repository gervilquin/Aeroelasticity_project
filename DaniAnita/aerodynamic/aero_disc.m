%% Aerodisc: Discretization for solving the aerodynamic problem
% INPUT
% y:     [1xlength(y)] the position in y of each STRUCTURAL node
% c:     chord of the profile
% xsc:   x position of the shear center
% OUTPUT:
% aerodynamic_points: [(length(y)-1)x3] position of each point in the
%                       aerodynamic center line
% collocation_points: [(length(y)-1)x3] position of each collocation point
% horse shoe:         Structure with [(length(y)-1)x3] matrices in fields
%                     A, B, C, and D. Contaning the position of the points
%                     defining the horseshoe associated to each aerodynamic
%                     point


function [collocation_points, aerodynamic_points,horse_shoe] = aero_disc(y,c,xsc,toplot)
    alpha = zeros(size(y));
    delta = zeros(size(y)); 
    
    for i=1:length(y)-1
        aerodynamic_points(i,1) = 0.25*c*cos(alpha(i));
        aerodynamic_points(i,2) = (y(i+1)+y(i))*0.5;
        aerodynamic_points(i,3) = delta(i) + 0.5*(delta(i+1)-delta(i))+...
            (xsc-0.25*c)*sin(alpha(i)+0.5*(alpha(i+1)-alpha(i)));
    end
    
    for i=1:length(y)-1
        collocation_points(i,1) = 0.75*c*cos(alpha(i));
        collocation_points(i,2) = (y(i+1)+y(i))*0.5;
        collocation_points(i,3) = delta(i) + 0.5*(delta(i+1)-delta(i))-...
            (0.75*c-xsc)*sin(alpha(i)+0.5*(alpha(i+1)-alpha(i)));
    end
    for i=1:length(y)-1

        horse_shoe.A(i,1) = 20*c*cos(alpha(i));
        horse_shoe.A(i,2) = y(i);
        horse_shoe.A(i,3) = delta(i)-(20*c-xsc)*sin(alpha(i));
        
        horse_shoe.B(i,1) = 0.25*c*cos(alpha(i));
        horse_shoe.B(i,2) = y(i);
        horse_shoe.B(i,3) = delta(i)+(xsc-0.25*c)*sin(alpha(i));
               

        horse_shoe.C(i,1) = 0.25*c*cos(alpha(i+1));
        horse_shoe.C(i,2) = y(i+1);
        horse_shoe.C(i,3) = delta(i+1)+(xsc-0.25*c)*sin(alpha(i+1));
        
        horse_shoe.D(i,1) = 20*c*cos(alpha(i+1));
        horse_shoe.D(i,2) = y(i+1);
        horse_shoe.D(i,3) = delta(i+1)-(20*c-xsc)*sin(alpha(i+1));
    end
    
    if toplot, plot_aero_disc(aerodynamic_points, collocation_points,...
            horse_shoe,xsc,y,delta); end
end
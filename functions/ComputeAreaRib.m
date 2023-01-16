function A_rib = ComputeAreaRib(t,c)
    A_rib = 0;

    % Computation without slot
    [x,airfoil_points] = CoordinatesSymmetricAirfoil(t,0.0075,0.095,100,c,'TE'); 
    for i=1:length(x)-1
        b = x(i+1)-x(i); h = airfoil_points(i+1)+airfoil_points(i);
        A = b*h;
        A_rib = A_rib + A;
    end
    
    % Substraction of slot
    b = 0.03; h = 0.005; A_slot = b*h;
    A_rib = A_rib - A_slot;
end
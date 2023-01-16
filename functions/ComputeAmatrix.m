function A = ComputeAmatrix(col_coor,panel_coor)
    
    % Initialize aerodynamic influece matrix
    A = zeros(size(col_coor,1)); 

    for i = 1:length(col_coor)
        for j = 1:length(col_coor)
            % Computation of distance vectors
            r1 = col_coor(i,:) - panel_coor(1,:,j);
            r2 = col_coor(i,:) - panel_coor(2,:,j);
            r3 = col_coor(i,:) - panel_coor(3,:,j);
            r4 = col_coor(i,:) - panel_coor(4,:,j);
            
            % Computation of vortex segment vectors
            l12 = panel_coor(2,:,j)-panel_coor(1,:,j);
            l23 = panel_coor(3,:,j)-panel_coor(2,:,j);
            l34 = panel_coor(4,:,j)-panel_coor(3,:,j);
            
            % Computation of induced velocities at collocation point
            v12 = (1/(4*pi)) * (cross(r1,r2))/(norm(cross(r1,r2))^2) * ( (dot(l12,r1)/norm(r1)) - (dot(l12,r2)/norm(r2)));
            v23 = (1/(4*pi)) * (cross(r2,r3))/(norm(cross(r2,r3))^2) * ( (dot(l23,r2)/norm(r2)) - (dot(l23,r3)/norm(r3)));
            v34 = (1/(4*pi)) * (cross(r3,r4))/(norm(cross(r3,r4))^2) * ( (dot(l34,r3)/norm(r3)) - (dot(l34,r4)/norm(r4)));
            
            % Computation of influence coefficient
            A(i,j) =dot((v12 + v23 + v34),[0,0,1]);
        end
    end
    

end
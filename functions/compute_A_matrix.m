function A = compute_A_matrix(col_coor,panel_coor)
    
    % Initialize the influece matrix
    A = zeros(length(col_coor));

    for i = 1:length(col_coor)
        for j = 1:length(col_coor)
            
            r_1 = col_coor(i,:) - panel_coor(1,:,j);
            r_2 = col_coor(i,:) - panel_coor(2,:,j);
            r_3 = col_coor(i,:) - panel_coor(3,:,j);
            r_4 = col_coor(i,:) - panel_coor(4,:,j);

%             l_12 = norm(panel_coor(1,:,j)-panel_coor(2,:,j));
%             l_23 = norm(panel_coor(2,:,j)-panel_coor(3,:,j));
%             l_34 = norm(panel_coor(3,:,j)-panel_coor(4,:,j));

            l_12 = panel_coor(2,:,j)-panel_coor(1,:,j);
            l_23 = panel_coor(3,:,j)-panel_coor(2,:,j);
            l_34 = panel_coor(4,:,j)-panel_coor(3,:,j);

            v_12 = (1/(4*pi)) * (cross(r_1,r_2))/(norm(cross(r_1,r_2))^2) * ( (dot(l_12,r_1)/norm(r_1)) - (dot(l_12,r_2)/norm(r_2)));
            v_23 = (1/(4*pi)) * (cross(r_2,r_3))/(norm(cross(r_2,r_3))^2) * ( (dot(l_23,r_2)/norm(r_2)) - (dot(l_23,r_3)/norm(r_3)));
            v_34 = (1/(4*pi)) * (cross(r_3,r_4))/(norm(cross(r_3,r_4))^2) * ( (dot(l_34,r_3)/norm(r_3)) - (dot(l_34,r_4)/norm(r_4)));

            
            A(i,j) =dot((v_12 + v_23 + v_34),[0,0,1]);
        end
    end
    

end
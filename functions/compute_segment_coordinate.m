function segment_coor = compute_segment_coordinate(y_nodal,ac_pos,chord_y)
    
    segment_coor = zeros(4,3,length(y_nodal)-1);

    for i = 1:length(y_nodal)-1
        
        segment_coor(1,:,i) = [20*chord_y(i),y_nodal(i),0];

        segment_coor(2,:,i) = [ac_pos(i,1),y_nodal(i),0];

        segment_coor(3,:,i) = [ac_pos(i,1),y_nodal(i+1),0];

        segment_coor(4,:,i) = [20*chord_y(i),y_nodal(i+1),0];

    end

end
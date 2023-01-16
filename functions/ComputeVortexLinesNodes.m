function segment_coor = ComputeVortexLinesNodes(y_el,ac_pos,c)
    %Initialisation
    nel = length(y_el)-1;
    segment_coor = zeros(4,3,nel); % [#vortex lines, ndim, nel]

    % Computation
    for i = 1:nel
        segment_coor(1,:,i) = [20*c,    y_el(i),    0];
        segment_coor(2,:,i) = [ac_pos,  y_el(i),    0];
        segment_coor(3,:,i) = [ac_pos,  y_el(i+1),  0];
        segment_coor(4,:,i) = [20*c,    y_el(i+1),  0];
    end
end
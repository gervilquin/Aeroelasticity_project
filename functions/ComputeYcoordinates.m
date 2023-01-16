function y_nodes = ComputeYcoordinates(y_sec,nsec,nesec)
% Compute nodal coordinates of the mesh elements.
    y_nodes = zeros((nsec*nesec)+1,1);
    for j = 1:nsec
        DeltaY = (y_sec(j+1)-y_sec(j))/nesec;
        for i = 1:nesec
            y_nodes(nesec*(j-1)+(i)) = y_sec(j) + DeltaY*(i-1);
        end
    end 
    y_nodes(end) = y_sec(end);
end


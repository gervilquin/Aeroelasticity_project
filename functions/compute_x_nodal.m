function x_nodal = compute_x_nodal(x_section,n_elements)
    % Compute the nodal coordinates between secitons 

    n_sections = length(x_section)-1;
    %x_nodal = zeros(1,n_element);
    x_nodal = [];

    for j = 1:n_sections

        delta_x = (x_section(j+1)-x_section(j))/n_elements;
        
        for i = 1:(n_elements)
            x_nodal(end+1) = x_section(j) + delta_x*(i-1);
        end

    end 
    x_nodal(end+1) = x_section(end);
end


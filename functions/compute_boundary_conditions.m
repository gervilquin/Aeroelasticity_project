function [Ip,If,u] = compute_boundary_conditions(Up,Ndof)
    
    u = zeros(Ndof,1);
    Ip = zeros(length(Up(:,1)),1);
    for i=1:length(Ip)
        Ip(i) = 6*(Up(i,2)-1) + Up(i,3);
        u(Ip(i),1) = Up(i,1);
    end
    If = setdiff(1:Ndof,Ip);

end
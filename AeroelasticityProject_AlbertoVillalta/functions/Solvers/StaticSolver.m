function [u] = StaticSolver(K,F,u,If,Ip)
    % System solver
    u(If,1) = K(If,If)\(F(If,1)-K(If,Ip)*u(Ip,1));
end
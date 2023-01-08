function solve_static_case(Nnod,y_nodal,I_au_0,I_au_1,I_au_2,I_fL,S,A_aero,K,U_inf,AoA_deg)

    Up = [  0 1 1;
            0 1 2;
            0 1 3]; % Prescribed DOF
    Ndof = Nnod*3;
    u_static = zeros(Ndof,1);
    Ip = zeros(length(Up(:,1)),1);
    for i=1:length(Ip)
        Ip(i) = 6*(Up(i,2)-1) + Up(i,3);
        u_static(Ip(i),1) = Up(i,1);
    end
    If = setdiff(1:Ndof,Ip);
    
    % compute force
    AoA = deg2rad(AoA_deg);
    u = AoA*ones(3*Nnod,1);
    u_dot = zeros(3*Nnod,1);
    u_dotdot = zeros(3*Nnod,1);
    
    alpha = I_au_0*u + I_au_1*u_dot + I_au_2*u_dotdot;
    L = -U_inf^2*S\A_aero*alpha;
    F= I_fL*L;
    
    % solve
    u_static(If) = K(If,If)\(F(If,1)-K(If,Ip)*u_static(Ip,1));
    
    % plot
    twist = u_static(1:3:Ndof-2);
    u_vertical = u_static(2:3:Ndof-1); 
    flection = u_static(3:3:Ndof);
    
    figure()
    subplot(3,1,1)
    plot(y_nodal,rad2deg(twist))
    ylabel("Twist","Interpreter","latex")
    subplot(3,1,2)
    plot(y_nodal,u_vertical)
    ylabel("U vertical","Interpreter","latex")
    subplot(3,1,3)
    plot(y_nodal,flection)
    ylabel("Deflection","Interpreter","latex")
    xlabel("Y [mm]",'Interpreter','latex')

end
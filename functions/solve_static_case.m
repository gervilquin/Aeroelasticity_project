function [twist,u_vertical,flection] = solve_static_case(Nnod,y_nodal,u_static,If,Ip,I_fL,S,A_aero,K,U_inf,rho_inf,AoA_deg,plot_bool)
    
%     Up = [  0 1 1;
%             0 1 2;
%             0 1 3]; % Prescribed DOF
%     
%     u_static = zeros(Ndof,1);
%     Ip = zeros(length(Up(:,1)),1);
%     for i=1:length(Ip)
%         Ip(i) = 6*(Up(i,2)-1) + Up(i,3);
%         u_static(Ip(i),1) = Up(i,1);
%     end
%     If = setdiff(1:Ndof,Ip);


    Ndof = Nnod*3;
    % compute force
    AoA = deg2rad(AoA_deg);
%     u = AoA*ones(3*Nnod,1);
%     u_dot = zeros(3*Nnod,1);
%     u_dotdot = zeros(3*Nnod,1);
    
    %alpha = I_au_0*u(If) + I_au_1*u_dot(If)  + I_au_2*u_dotdot(If) ;
    alpha = AoA*ones(length(If)/3,1);
    S_aero = -U_inf^2*rho_inf*S;
    L = S_aero*inv(A_aero)*alpha;
    F= I_fL*L;
    
    % solve
    u_static(If) = K(If,If)\(F-K(If,Ip)*u_static(Ip,1));
    
    % plot
    twist = u_static(1:3:Ndof-2);
    u_vertical = u_static(2:3:Ndof-1); 
    flection = u_static(3:3:Ndof);
    
    if plot_bool == true
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

end
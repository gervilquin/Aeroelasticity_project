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
        ylabel("Twist [deg]","Interpreter","latex")
        title(strcat("U_{\infty} = ",num2str(U_inf)," m/s   AoA = ",num2str(AoA_deg)," deg"))
        grid on
        grid minor
        subplot(3,1,2)
        plot(y_nodal,u_vertical)
        ylabel("U vertical [m]","Interpreter","latex")       
        grid on
        grid minor
        subplot(3,1,3)
        plot(y_nodal,flection)
        ylabel("Deflection","Interpreter","latex")
        xlabel("Y [mm]",'Interpreter','latex')
        grid on
        grid minor

        % Theoretical lift
        L_t = pi/(1+2*pi/((2*0.55)^2/(0.1*2*0.55)))*AoA*rho_inf*0.55*0.1*U_inf^2;
        disp(['Computed with A matrix = ',num2str(sum(L))])
        disp(['Theroetical lift       = ',num2str(L_t)])
        disp(['Lift after I_fl = ',num2str(sum(F(2:3:end)))])
        
        figure()
        hold on
        plot(y_nodal(1:end-1),L)
        plot(y_nodal(1:end-1),L_t/length(L)*ones(length(L)))
        xlabel("Y [m]","Interpreter","latex")
        ylabel("Lift [N/m]")
        legend(["Numerical Lifting line","Theoretical lifting line"])
        hold off
    end

end
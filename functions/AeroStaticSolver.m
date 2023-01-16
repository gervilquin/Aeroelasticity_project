function [twist,u_vertical,flection] = AeroStaticSolver(y_el,u,If,Ip,I_fL,I_au_0,S,A_aero,K,U_inf,rho_inf,AoA_deg,plot_bool)
    % Initialisation 
    nel = length(y_el)-1;
    max_tol = 1e-6;
    max_error = 1;
    u_ite = zeros(length(u),1);
    ite = 1;

    % Compute force vector
    AoA = deg2rad(AoA_deg);
    alpha = AoA*ones(nel,1);
    S_aero = -(U_inf^2)*rho_inf*S;
%     L = S_aero*inv(A_aero)*alpha;
%     F= I_fL*L;

    % Start loop to find the equilibrium of the system.
    while (ite <= 100) && (max_error >= max_tol)
        
        % Compute the lift force
        L = S_aero*inv(A_aero)*alpha;
        F= I_fL*L;

        % System solver
        u(If,1) = K(If,If)\(F(If,1)-K(If,Ip)*u(Ip,1));
        
        % Compute max error
        max_error = max(abs(u-u_ite));

        %if max_error > max_tol
        alpha = I_au_0*u;
        %end

        u_ite = u;

        ite = ite + 1;
    end


    % Results presentation
    twist = u(1:3:end);
    u_vertical = u(2:3:end); 
    flection = u(3:3:end);
    
    if plot_bool == true
        figure()
        subplot(3,1,1)
        hold on
        plot(y_el,rad2deg(twist))
        plot(y_el(2:end),rad2deg(alpha))
        ylabel("Twist [deg]","Interpreter","latex")
        title(strcat("U_{\infty} = ",num2str(U_inf)," m/s   AoA = ",num2str(AoA_deg)," deg"))
        grid on
        grid minor
        legend(["twist","alpha"])
        hold off
        subplot(3,1,2)
        plot(y_el,u_vertical)
        ylabel("U vertical [m]","Interpreter","latex")       
        grid on
        grid minor
        subplot(3,1,3)
        plot(y_el,flection)
        ylabel("Deflection","Interpreter","latex")
        xlabel("Y [mm]",'Interpreter','latex')
        grid on
        grid minor

        % Theoretical lift
%         L_t = pi/(1+2*pi/((2*0.55)^2/(0.1*2*0.55)))*AoA*rho_inf*0.55*0.1*U_inf^2;
%         disp(['Computed with A matrix = ',num2str(sum(L))])
%         disp(['Theroetical lift       = ',num2str(L_t)])
%         disp(['Lift after I_fl = ',num2str(sum(F(2:3:end)))])
%         
%         figure()
%         hold on
%         plot(y_el(1:end-1),L)
%         plot(y_el(1:end-1),L_t/length(L)*ones(length(L)))
%         xlabel("Y [m]","Interpreter","latex")
%         ylabel("Lift [N/m]")
%         legend(["Numerical Lifting line","Theoretical lifting line"])
%         hold off
    end
end
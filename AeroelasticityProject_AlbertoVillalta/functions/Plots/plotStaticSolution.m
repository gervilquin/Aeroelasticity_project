function plotStaticSolution(u,F,L,y_el,U_inf,rho_inf,AoA_deg,plot_bool)

   % Resultat displacements for each DOF
   twist = u(1:3:end);
   u_vertical = u(2:3:end); 
   deflection = u(3:3:end);
    
   % Results presentation
   if plot_bool == true
        figure(1)
        subplot(3,1,1)
        plot(y_el,rad2deg(twist))

        ylim([AoA_deg,AoA_deg+1.1*(max(rad2deg(twist))-AoA_deg)])

        ylabel("$\theta$ [deg]","Interpreter","latex")
        title(strcat("$U_{\infty}$ = ",num2str(U_inf)," m/s   AoA = ",num2str(AoA_deg)," deg"),'Interpreter','latex')
        grid on
        grid minor
        subplot(3,1,2)
        plot(y_el,u_vertical)
        ylim([0,1.1*max(u_vertical)])
        ylabel("$w_{sc}$[m]","Interpreter","latex")       
        grid on
        grid minor
        subplot(3,1,3)
        plot(y_el,deflection)
        ylim([0,1.1*max(deflection)])
        ylabel("$\gamma$ [-]","Interpreter","latex")
        xlabel("Y [m]",'Interpreter','latex')
        grid on
        grid minor

        % Theoretical lift
        L_t = pi/(1+2*pi/((2*0.55)^2/(0.1*2*0.55)))*deg2rad(AoA_deg)*rho_inf*0.55*0.1*U_inf^2;
        disp(['Computed with A matrix = ',num2str(sum(L))])
        disp(['Theroetical lift       = ',num2str(L_t)])
        disp(['Lift after I_fl = ',num2str(sum(F(2:3:end)))])
        
%         figure(2)
%         hold on
%         plot(y_el(1:end-1),L)
%         plot(y_el(1:end-1),L_t/length(L)*ones(length(L)))
%         xlabel("Y [m]","Interpreter","latex")
%         ylabel("Lift [N/m]")
%         legend(["Numerical Lifting line","Theoretical lifting line"])
%         hold off
    end
end
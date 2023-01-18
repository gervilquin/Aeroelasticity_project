function plot2StaticSolution(u,F,L,y_el,U_inf,rho_inf,AoA_deg,plot_bool)

   % Resultat displacements for each DOF
   twist1 = u(:,1);
   u_vertical1 = u(:,2); 
   deflection1 = u(:,3);

   twist2 = u(:,4);
   u_vertical2 = u(:,5); 
   deflection2 = u(:,6);
    
   % Results presentation
   if plot_bool == true
        figure()
        subplot(3,1,1)
        hold on
        plot(y_el,rad2deg(twist1))
        plot(y_el,rad2deg(twist2))
        ylim([AoA_deg,AoA_deg+1.1*(max(rad2deg(twist2))-AoA_deg)])
        ylabel("$\theta$ [deg]","Interpreter","latex")
        title(strcat("$U_{\infty}$ = ",num2str(U_inf)," m/s   AoA = ",num2str(AoA_deg)," deg"),'Interpreter','latex')
        grid on
        grid minor
        hold off

        subplot(3,1,2)
        hold on
        plot(y_el,u_vertical1)
        plot(y_el,u_vertical2)
        ylim([0,1.1*max(u_vertical2)])
        legend(["Static","Aerodynamic coupling"],'interpreter','latex','Location','northwest')
        ylabel("$w_{sc}$[m]","Interpreter","latex")       
        grid on
        grid minor
        hold off


        subplot(3,1,3)
        hold on
        plot(y_el,deflection1)
        plot(y_el,deflection2)
        ylim([0,1.1*max(deflection2)])
        ylabel("$\gamma$ [-]","Interpreter","latex")
        xlabel("Y [m]",'Interpreter','latex')
        grid on
        grid minor
        hold off

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
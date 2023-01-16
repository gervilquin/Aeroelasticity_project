function plot_vars_static(aerodynamic_point,y,u,f,alpa,L)
            figure(1)
            subplot(2,4,1);
            hold on
            grid on
            plot(aerodynamic_point(:,2),L,'LineWidth',1.5)
            xlabel('y')
            ylabel('L')
    
            subplot(2,4,2)
            hold on
            grid on 
            ylabel('Force over vertical deflection \omega_i')
            xlabel('y')
            plot(y,f(2:3:end),'LineWidth',1.5)
            
            subplot(2,4,3)
            hold on
            grid on
            ylabel('Moment over torsion angle \Theta_i')
            xlabel('y')
            plot(y,f(1:3:end),'LineWidth',1.5)
    
            subplot(2,4,4)
            hold on
            grid on
            ylabel('Moment over bending angle \gamma_i')
            xlabel('y')
            plot(y,f(3:3:end),'LineWidth',1.5)

            subplot(2,4,6)
            hold on
            grid on
            ylabel('\omega_i')
            xlabel('y')
            plot(y,u(2:3:end),'LineWidth',1.5)

            subplot(2,4,7)
            hold on
            grid on
            ylabel('\Theta_i')
            xlabel('y')
            plot(y,rad2deg(u(1:3:end)),'LineWidth',1.5)

            subplot(2,4,8)
            hold on
            grid on
            ylabel('\gamma_i')
            xlabel('y')
            plot(y,rad2deg(u(3:3:end)),'LineWidth',1.5)

            subplot(2,4,5)
            hold on
            grid on
            ylabel('\alpha^*')
            xlabel('y')
            plot(aerodynamic_point(:,2),rad2deg(alpa),'LineWidth',1.5)

end
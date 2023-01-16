function  plot_deformations(y,u)
    figure
    hold on
    grid minor
    ylabel('Vertical deflection \omega_i')
    xlabel('y')
    plot(y,u(2:3:end),'LineWidth',1.5)
    
    figure
    hold on
    grid minor
    ylabel('Torsion angle \Theta_i')
    xlabel('y')
    plot(y,rad2deg(u(1:3:end)),'LineWidth',1.5)
    
    figure
    hold on
    grid minor
    ylabel('Bending angle \gamma_i')
    xlabel('y')
    plot(y,rad2deg(u(3:3:end)),'LineWidth',1.5)
end
function  plot_modal(y,phi,w2,modes)
    phi = real(phi);   
    
    
    figure
    subplot(1,3,1)  
    hold on
    grid minor
    for i=1:modes
        plot(1e3*y,rad2deg(phi(1:3:end,i))./max(rad2deg(phi(1:3:end,i))),'LineWidth',1.5)
        leg{i} = strcat('f=',num2str(sqrt(w2(i))/2/pi));
    end
    xlabel('y [mm]')
    ylabel('\Theta [deg]')
    legend(leg)

    subplot(1,3,2)  
    hold on
    grid minor

    for i=1:modes
        plot(1e3*y,1e3*(phi(2:3:end,i))./max(1e3*(phi(2:3:end,i))),'LineWidth',1.5)
    end
    xlabel('y [mm]')
    ylabel('\omega [mm]')
    legend(leg)

    subplot(1,3,3)  
    hold on
    grid minor

    for i=1:modes
        plot(1e3*y,rad2deg(phi(3:3:end,i))./max(rad2deg(phi(3:3:end,i))),'LineWidth',1.5)
    end
    xlabel('y [mm]')
    ylabel('\gamma [deg]')
    legend(leg)
end
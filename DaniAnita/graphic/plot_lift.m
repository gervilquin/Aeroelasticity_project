function plot_lift(L,collocation_points)
    figure
    hold on
    grid on
    grid minor
    plot(1e3.*collocation_points(:,2),L,'o-','LineWidth',1.5);
    title('Lift along y direction')
    ylabel('L [N]')
    xlabel('y [mm]');
end
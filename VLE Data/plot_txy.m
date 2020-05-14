function plot_txy()
    load('base_xy_diagram.mat','xEtoh','yEtoh','T_eqlm')
    
    plot(xEtoh,T_eqlm,'b','LineWidth',2)
    hold on
    plot(yEtoh,T_eqlm,'g','LineWidth',2)
    xlabel('xEtoh,yEtoh')
    ylabel(['Temperature ' char(176)  'C'])
    saveas(gcf,'xy diagrams/T_xy_1atm.png')
    close(gcf)
end
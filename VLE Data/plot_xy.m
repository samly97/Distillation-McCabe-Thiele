function plot_xy()
    load('base_xy_diagram.mat','x','y','xEtoh','yEtoh')
    
    plot(xEtoh,yEtoh,'b','LineWidth',2)
    hold on
    plot(x,y,'g','LineWidth',2)
    xlabel('x, Mole Fraction Ethanol')
    ylabel('y, Mole Fraction Ethanol')
    saveas(gcf,'xy diagrams/base_xy_1atm.png')
    close(gcf)
end
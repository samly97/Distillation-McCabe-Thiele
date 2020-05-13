function plot_xy(n)
    vle = generate_vle_data(n);
    xEtoh=vle(:,1); yEtoh=vle(:,2); T_eqlm =vle(:,3);
    x=linspace(0,1,n)'; y=x;
    save('base_xy_diagram.mat','x','y','xEtoh','yEtoh','T_eqlm')
    
    plot(xEtoh,yEtoh,'b','LineWidth',2)
    hold on
    plot(x,y,'g','LineWidth',2)
    xlabel('x, Mole Fraction Ethanol')
    ylabel('y, Mole Fraction Ethanol')
    saveas(gcf,'base_xy_1atm.png')
end
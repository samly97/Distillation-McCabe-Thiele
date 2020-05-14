% input:
%   Ae: Area fraction of ethanol from GC - in order from feed, distillate,
%   and bottom
%   T: Temperature of the feed
%   R: Reflux ratio
%   F: Feed flowrate (mol/time)
%   D: Distillate flowrate (mol/time)
function main(Ae,T,R,F,D)
    % Find paths
    Project_Folder = pwd;
    VLE = sprintf('%s/VLE Data',Project_Folder);

    % Add dependencies
    addpath(genpath(VLE))
    
    % find ethanol-water equilibrium data
    load_vle_data(200)
    
    % plot xy and Txy diagrams
    plot_xy()
    plot_txy()
    
    % McCabe Thiele
    temp = gc2frac(Ae);
    x = temp(:,2)';
    temp = find_q(x(1),T);
    q = temp(3);
    lewis_sorel(x, R, q, F, D)
    
    % remove files created in load_vle_data()
    cleanup()
end

function load_vle_data(n)
    % find liquid-vapour equilibrium ethanol concentration
    vle = generate_vle_data(n);
    xEtoh=vle(:,1); yEtoh=vle(:,2); T_eqlm =vle(:,3);
        x=linspace(0,1,n)'; y=x;
    save('base_xy_diagram.mat','x','y','xEtoh','yEtoh','T_eqlm')
    eqlm_boiling = [vle(:,1) vle(:,3)];
    save('eqlm_boiling.mat','eqlm_boiling')
end

function cleanup()
    delete('base_xy_diagram.mat')
    delete('eqlm_boiling.mat')
end
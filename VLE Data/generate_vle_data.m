% inputs:
%   n: how many eqlm points to generate

% output:
%   n x 3 matrix with x_eqlm, y_eqlm, and T_eqlm as column headers
function xy_eqlm = generate_vle_data(n)
    T_bp_e = 79.292255;
    T_bp_w = 101.981615;
    C = setup();
    xy_eqlm = find_y_eqlm(C,[T_bp_e T_bp_w],n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating VLE Data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = find_y_eqlm(C,T_bounds,n)
    data = setup_eqlm_data(T_bounds,n);
    for i = 2:n-1
        x = data(i,1); To = data(i-1,3);
        T = find_equilibrium(C{1},x,To);
        data(i,3) = T;
        data(i,2) = C{2}(x,T);
    end
end

function T_eqlm = find_equilibrium(P,x,To)
    P_tar = @(T) P(x,T)-760;
    T_eqlm = fzero(P_tar,To);
end

function data = setup_eqlm_data(T_bounds,n)
    x = linspace(0,1,n)';
    data = zeros(n,3);
    data(:,1) = x;
    data(1,3)=T_bounds(2); % pure ethanol
    data(n,2)=1; data(n,3)=T_bounds(1); % pure water
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup Functions Which Relate Liquid Composition to Total Pressure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = setup()
    ant_e = [8.1348 1662.484 238.1311];
    ant_w = [8.0557 1723.643 233.0764];
    G12 = 0.1388;
    G21 = 0.9296;

    P_v = @(T, A, B, C) 10.^(A-B./(C+T));
    P_ve = @(T) P_v(T, ant_e(1), ant_e(2), ant_e(3));
    P_vw = @(T) P_v(T, ant_w(1), ant_w(2), ant_w(3));
    x2 = @(x1) 1-x1;
    ln_gamma1 = @(x1) -log(x1+x2(x1).*G12)+x2(x1)*(G12./(x1+x2(x1).*G12)-G21./(x2(x1)+x1.*G21));
    ln_gamma2 = @(x1) -log(x2(x1)+x1.*G21)-x1*(G12./(x1+x2(x1).*G12)-G21./(x2(x1)+x1.*G21));
    
    P = @(x1,T) x1.*exp(ln_gamma1(x1)).*P_ve(T) + x2(x1).*exp(ln_gamma2(x1)).*P_vw(T);
    y1 = @(x1,T) x1.*exp(ln_gamma1(x1)).*P_ve(T)/760;
    C = {P y1};
end
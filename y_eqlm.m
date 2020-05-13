% inputs:
%   xi: as a n x 1 vector

% assumptions:
%   fzero will solve based on in-between point of boiling temperatures

function y_eqlm = y_eqlm(xi)
    T_bp_e = 79.292255;
    T_bp_w = 101.981615;
    To = (T_bp_e + T_bp_w)/2;
    C = setup();
    T_eqlm = find_equilibrium(C{1},xi,To);
    y_eqlm = C{2}(xi,T_eqlm);
end

function T_eqlm = find_equilibrium(P,x,To)
    T_eqlm = zeros(length(x),1);
    for i = 1:length(x)
        P_tar = @(T) P(x(i),T)-760;
        T_eqlm(i) = fzero(P_tar,To);
    end
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
    ln_gamma1 = @(x1) -log(x1+x2(x1).*G12)+x2(x1).*(G12./(x1+x2(x1).*G12)-G21./(x2(x1)+x1.*G21));
    ln_gamma2 = @(x1) -log(x2(x1)+x1.*G21)-x1.*(G12./(x1+x2(x1).*G12)-G21./(x2(x1)+x1.*G21));
    
    P = @(x1,T) x1.*exp(ln_gamma1(x1)).*P_ve(T) + x2(x1).*exp(ln_gamma2(x1)).*P_vw(T);
    y1 = @(x1,T) x1.*exp(ln_gamma1(x1)).*P_ve(T)/760;
    C = {P y1};
end
% inputs:
%   Ae: matrix or scalar, area of ethanol
%
% output:
%   n x 2 matrix, 1st column mass frac and 2nd mol frac

function fractions = gc2frac(Ae)
    C = setup();
    fractions = calc_fracs(C,Ae);
end

function data = calc_fracs(C,Ae)
    n = length(Ae);
    data = zeros(n,2);
    
    for i = 1:n
        data(i,1)=C{1}(Ae(i));
        data(i,2)=C{2}(Ae(i));
    end
end

function C = setup()
    w_e=0.64; w_w=0.55; Me=46; Mw=18;

    Aw = @(Ae) 1-Ae;
    xm_e = @(Ae) Ae.*w_e/(Ae.*w_e+Aw(Ae).*w_w);
    xm_w = @(Ae) 1-xm_e(Ae);
    x_e = @(Ae) xm_e(Ae)./Me./(xm_e(Ae)./Me+xm_w(Ae)/Mw);
    C = {xm_e x_e};
end
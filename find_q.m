% scripts assumes scalar 
%
% input:
%   x: mol fraction ethanol in feed
%   T: feed temperature (celsius)
%
% output:
%   xe, Tf, q, Tb, Sens_avg, Sens_e, Sens_w, Hv_avg, Hv_e, Hv_w
function data = find_q(x,Tf)
    C = setup();
    data = generate_data(C,x,Tf);
end

%%% Calculate individual quantities then store
function data = generate_data(C,x,Tf)
    Tf = Tf + 273;
    Tb = C{1}(x);

    data = zeros(1,10);
    data(1,1)=x; data(1,2)=Tf-273; data(1,4)=Tb-273; % x, Tf, Tb
    data(1,6)=C{2}(Tf,Tb); data(1,7)=C{3}(Tf,Tb); % sens_e, sens_w
    data(1,9)=C{4}(Tb); data(1,10)=C{5}(Tb); % Hv_e, Hv_w
    data(1,5)=x*data(1,6)+(1-x)*data(1,7); % sens_avg
    data(1,8)=x*data(1,9)+(1-x)*data(1,10); % Hv_avg
    data(1,3)=(data(1,8)+data(1,5))/data(1,8); % q
end

%%%%%%%%%%%%%
%%% Setup %%%
%%%%%%%%%%%%%
function C = setup()
    Tx_eqlm = setup_Txy_eqlm(); % K
    Sens = setup_sens(); % J/mol
    Hv = setup_latent(); % J/mol
    C = {Tx_eqlm Sens{1} Sens{2} Hv{1} Hv{2}};
end

function Txy_eqlm = setup_Txy_eqlm()
    %%% Txy Equilibrium %%%
    load eqlm_boiling.mat
    xe=eqlm_boiling(:,1)'; Tb=eqlm_boiling(:,2)';
    b = lin_reg(xe,Tb,8);
    Txy_eqlm = @(x) b(1) + b(2).*x + b(3).*x.^2 + b(4).*x.^3 + b(5).*x.^4+ ...
    b(6).*x.^5 + b(7).*x.^6 + b(8).*x.^7 + b(9).*x.^8; % C
    Txy_eqlm = @(x) Txy_eqlm(x) + 273; % K % final eqn
end

function sens = setup_sens()
    %%% Sensible Heat %%%
    % From Perry's 1999
    Cp = [102640 -139.63  -0.030341 0.0020386      0            % ethanol
          276370 -2090.10 8.125     -1.4116*10^-2  9.37*10^-6]; % water
    % integrate from T1 to T2
    Cp_e = @(T) Cp(1,1)+Cp(1,2).*T+Cp(1,3).*T.^2+Cp(1,4).*T.^3+Cp(1,5).*T.^4;
    Cp_w = @(T) Cp(2,1)+Cp(2,2).*T+Cp(2,3).*T.^2+Cp(2,4).*T.^3+Cp(2,5).*T.^4;
    Sens_Cp_e = @(T1,T2) integral(Cp_e,T1,T2)/1000; % J/mol
    Sens_Cp_w = @(T1,T2) integral(Cp_w,T1,T2)/1000; % J/mol
    sens = {Sens_Cp_e Sens_Cp_w};
end

function Hv_fun = setup_latent()
    % T is feed temperature
    %%% Latent Heat %%%
    A = [60.8036 54]; % eth wat
    Tc = [516.25 647.13]; % eth wat
    n = [0.38 0.34]; % eth wat
    Hv = @(T,A,Tc,n) A*(1-T./Tc).^n*1000;
    Hv_e = @(T) Hv(T,A(1),Tc(1),n(1)); % J/mol % final eqn
    Hv_w = @(T) Hv(T,A(2),Tc(2),n(2)); % J/mol % final eqn
    Hv_fun = {Hv_e Hv_w};
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% linear regression %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
function out = lin_reg(x, y, n) 
    for i = 1:n+1
        X(:,i) = x.^(i-1)';
    end
    A = X'*X;
    b = X'*y';
    out = A\b;
end
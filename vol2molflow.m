% inputs:
%   V_flow: Feed/Bottom/Distillate flowrate in L/h (n x 1)
%   xe: Feed/Bottom/Distillate mol frac ethanol (n x 1)
% output:
%   mol_flow: Given flowrate in mol/h (n x 1)

% assumptions:
%   T and P at 25 C and 1 atm, respectively
function mol_flow = vol2molflow(V_flow,xe)    
    f = setup;
    V_soln = sol_volume(f,xe);
    mol_flow = convert2mol(V_flow,V_soln);
end

%%%%%%%%%%%%%
%%% Setup %%%
%%%%%%%%%%%%%
function f = setup()
    Vmix_n = @(x) -46.283.*x.^6 + 137.55.*x.^5 - 147.62.*x.^4 + ...
        64.952.*x.^3 - 2.993.*x.^2 - 5.6514.*x + 0.0292;
    f = Vmix_n; % ml/mol
end

%%%%%%%%%%%%%%%%%%%%%%
%%% Molar Flowrate %%%
%%%%%%%%%%%%%%%%%%%%%%
function mol_flow = convert2mol(V_flow,V_soln)
    V_soln = V_soln/1000; % L/h
    mol_flow = V_flow./V_soln;
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% Solution Volume %%%
%%%%%%%%%%%%%%%%%%%%%%%
function V_soln = sol_volume(f,x)
    par_vol = partial_volume(f,x);
    V_soln = x.*par_vol(:,1) + (1-x).*par_vol(:,2); % ml/mol
end

%%%%%%%%%%%%%%%%%%%%%%
%%% Partial Volume %%%
%%%%%%%%%%%%%%%%%%%%%%
% out: V_eth V_wat
function out = partial_volume(f,x)
    Ve_pure = 1000/789*46.07; % ml/mol
    Vw_pure = 1000/997*18.02; % ml/mol
    
    % y = mx + b
    % int:0 --> Vw - Vw_pure
    % int:1 --> Ve - Ve_pure
    m = differentiate(f,x);
    b = f(x) - m.*x;
    
    Ve = (m+b) + Ve_pure;
    Vw = b + Vw_pure;
    out = [Ve Vw];
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%% First Derivative %%%
%%%%%%%%%%%%%%%%%%%%%%%%
% h: step
% First Difference w O(h^2) error
function deriv = differentiate(f,x)
    h = 0.001; 
    deriv = (-f(x+2.*h) + 4*f(x+h) - 3*f(x))./(2*h);
end
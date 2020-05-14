% inputs:
%   comp: specified composition in mole fraction ethanol: x_f, x_d, x_b
%   R: reflux ratio
%   q: feed quality

% assumptions:
%   constant molal overflow based on similar latent heat of vaporization

function lewis_sorel(comp, R_in, q_in, F_in, D_in)
    setup(comp, R_in, q_in, F_in, D_in)
    plot_distillation()
    theoretical_trays();
    
    fig_title = sprintf('VLE for Ethanol-Water at 1atm (R =%.2f)',R_in);
    xlabel('x, Mol Fraction Ethanol')
    ylabel('y, Mol Fraction Ethanol')
    title(fig_title)
    
    file_loc = sprintf('McCabe Thiele Diagrams/VLE Eth Wat 1atm R %.2f.png',R_in);
    saveas(gcf,file_loc)
    close(gcf)
end

%%%%%%%%%%%%%
%%% setup %%%
%%%%%%%%%%%%%
function setup(comp, R_in, q_in, F_in, D_in) 
    setup_distill_params(comp, R_in, q_in, F_in, D_in)
    setup_eqns()
end

function setup_distill_params(comp, R_in, q_in, F_in, D_in)
    global x_f x_d x_b q F D R L L_bar V_bar
    x_f = comp(1); x_d = comp(2); x_b = comp(3);
    q = q_in; F = F_in; D = D_in; R = R_in;
    L = R*D;
    temp = calc_strip_params(q,F,L);
    L_bar = temp(1); V_bar = temp(2);
end

function setup_eqns()
    global eqns q x_f L_bar V_bar x_b L D x_d
    load base_xy_diagram.mat
    b = lin_reg(xEtoh',yEtoh',6);
    q_b = q_line(q, x_f);  
    r_b = rect_line(L, D, x_d);
    s_b = strip_line();
    
    eqlm = @(x) b(1) + b(2)*x + b(3)*x.^2 + b(4)*x.^3 + b(5)*x.^4 + b(6)*x.^5 + b(7)*x.^6;
    q_fit = @(x) q_b(2) + q_b(1).*x;
    s_fit = @(x) s_b(2) + s_b(1).*x;
    r_fit = @(x) r_b(2) + r_b(1).*x;
    eqns = {};
    eqns = {eqlm, q_fit, s_fit, r_fit};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate stripping params %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = calc_strip_params(q,F,L)
    L_bar = q*F + L;
    V_bar = L-F+(L_bar - L);
    out = [L_bar, V_bar];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% feed, stripping, and rectifying line coefficients %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = q_line(q, x_f)
    m = q/(q-1);
    b = 1/(1-q)*x_f;
    out = [m, b];
end
function out = strip_line()
    global L D x_d q x_f x_b
    r_b = rect_line(L,D,x_d);
    q_b = q_line(q,x_f);
    f = @(x) (r_b(1).*x+r_b(2))-(q_b(1).*x+q_b(2));
    nice = fzero(f,x_f);
    r_line = @(x) r_b(1).*x+r_b(2);
    m = (r_line(nice)-x_b)/(nice-x_b);
    b = x_b - m*x_b;
    out = [m, b];
end

function out = rect_line(L, D, x_d)
    m = L/(L+D);
    b = D/(L+D)*x_d;
    out = [m, b];
end

%%%%%%%%%%%%%%%%
%%% plotting %%%
%%%%%%%%%%%%%%%%
function plot_distillation()
    global eqns x_f x_d x_b
    load base_xy_diagram.mat
    
    plot(xEtoh, yEtoh,'b',x,y,'k','LineWidth',3)
    hold on
    plot(x_b,x_b,'or','MarkerSize',10)
    hold on
    x_q = linspace(x_f,eql_intercept(x_f,'q'),16);
    plot(x_q,eqns{2}(x_q),'m','LineWidth',2)
    hold on 
    x_r = linspace(x_d,q_int(x_f),16);
    plot(x_r,eqns{4}(x_r),'g','LineWidth',2);
    hold on
    x_s = linspace(x_b,q_int(x_f),16);
    plot(x_s,eqns{3}(x_s),'c','LineWidth',2);
    hold on
end

function theoretical_trays()
    global eqns x_f x_d x_b
    draw_tray(eqns, x_d, x_f, x_b);
end

function draw_tray(eqns, x, x_f, x_b)
    x_eqlm = eql_intercept(x,'');
    y_eqlm = eqns{1}(x_eqlm);
    f = current_section(x);
    xp = linspace(x,x_eqlm,20);
    yp = linspace(f(x),f(x),20);
    plot(xp,yp,'-.b');
    hold on
    
    f = current_section(x_eqlm);
    xp = linspace(x_eqlm,x_eqlm,20);
    yp = linspace(y_eqlm,f(x_eqlm),20);
    plot(xp,yp,'-.b')
    hold on
    
    if x_eqlm > x_b
       draw_tray(eqns,x_eqlm,x_f,x_b);
    end
end

function f = current_section(x_op)
    global x_f eqns
    if x_op > x_f
        f = @(x) eqns{4}(x);
    else
        f = @(x) eqns{3}(x);
    end
end

function out = eql_intercept(xo,line)
    global eqns
    if line == 'q'
        f = @(x) eqns{1}(x) - eqns{2}(x);
    else
        e = current_section(xo);
        f = @(x) eqns{1}(x) - e(xo);
    end
    out = fzero(f,xo);
end

function out = q_int(xo)
    global eqns
    f = @(x) eqns{4}(x) - eqns{2}(x);
    out = fzero(f,xo);
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
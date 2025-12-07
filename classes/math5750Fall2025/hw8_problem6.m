%% Slow (v,n) subsystem of fast-slow Hodgkin-Huxley
% Phase plane, nullclines, trajectories, manifolds, and bifurcation vs I_app

clear; clc; close all;

%% Parameters
gNa = 120;
gK  = 36;
gL  = 0.3;
Cm  = 1;

vNa = 115;
vK  = -12;
vL  = 10.6;

% Quasi-steady Na activation m = m_inf(v)
alpha_m = @(v) alpha_m_fun(v);          % safe version below
beta_m  = @(v) 4*exp(-v/18);
m_inf   = @(v) alpha_m(v) ./ (alpha_m(v) + beta_m(v));

% Slow K gate n
alpha_n = @(v) alpha_n_fun(v);          % safe version below
beta_n  = @(v) 0.125*exp(-v/80);

% Rest-like initial condition
v0 = -65;
n0 = 0.3176;
y0 = [v0; n0];

%% 1. Phase plane, nullclines, and trajectory for one I_app
Iapp_pp = 6;   % adjust to explore how the phase plane changes

% v-n grid
vmin = -80; vmax = 120;  nv = 150;
nmin = 0;   nmax = 1;    nn = 150;
[v_grid, n_grid] = meshgrid(linspace(vmin, vmax, nv), ...
                            linspace(nmin, nmax, nn));

[vdot, ndot] = slow_vector_field(v_grid, n_grid, Iapp_pp, ...
                                 gNa,gK,gL,Cm,vNa,vK,vL, ...
                                 m_inf,alpha_n,beta_n);

figure(1); clf; hold on;
quiver(v_grid, n_grid, vdot, ndot, 0.7, 'k', ...
       'DisplayName','vector field');
xlabel('v (mV)');
ylabel('n');
title(sprintf('Slow (v,n) phase plane, I_{app} = %.2f', Iapp_pp));
box on;

% Nullclines: v-dot = 0 (red), n-dot = 0 (blue)
contour(v_grid, n_grid, vdot, [0 0], 'r', 'LineWidth', 2, ...
        'DisplayName','v-nullcline');
contour(v_grid, n_grid, ndot, [0 0], 'b', 'LineWidth', 2, ...
        'DisplayName','n-nullcline');

% Trajectory
tspan = [0 2000];
ode_fun = @(t,y) slow_rhs(t,y,Iapp_pp, ...
                          gNa,gK,gL,Cm,vNa,vK,vL, ...
                          m_inf,alpha_n,beta_n);
[ts, y_traj] = ode15s(ode_fun, tspan, y0);
plot(y_traj(:,1), y_traj(:,2), 'k', 'LineWidth', 2, ...
     'DisplayName','trajectory');

% Time series
figure(2); clf;
subplot(2,1,1);
plot(ts, y_traj(:,1), 'LineWidth', 1.5);
ylabel('v (mV)');
title(sprintf('Trajectory for I_{app} = %.2f', Iapp_pp));

subplot(2,1,2);
plot(ts, y_traj(:,2), 'LineWidth', 1.5);
xlabel('t (ms)');
ylabel('n');

%% 2. Equilibrium and manifolds at this I_app (only if saddle)
F_eq = @(Y) slow_rhs(0,Y,Iapp_pp, ...
                     gNa,gK,gL,Cm,vNa,vK,vL, ...
                     m_inf,alpha_n,beta_n);
Yguess = [-60; 0.3];

try
    opts = optimoptions('fsolve','Display','off');
    [Ystar,~,exitflag] = fsolve(F_eq, Yguess, opts);
catch
    exitflag = -1;
end

if exitflag > 0
    v_star = Ystar(1);
    n_star = Ystar(2);

    figure(1); hold on;
    plot(v_star, n_star, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8, ...
         'DisplayName','equilibrium');

    % Numerical Jacobian
    J = numerical_jacobian(@(Y) slow_rhs(0,Y,Iapp_pp, ...
                               gNa,gK,gL,Cm,vNa,vK,vL, ...
                               m_inf,alpha_n,beta_n), Ystar);
    [V_eig, D_eig] = eig(J);
    lambda = diag(D_eig);

    if all(imag(lambda) == 0) && prod(sign(real(lambda))) < 0
        eps_man = 1e-3;

        [~, i_unst] = max(real(lambda));
        [~, i_stab] = min(real(lambda));
        v_unst = V_eig(:, i_unst);
        v_stab = V_eig(:, i_stab);

        % Unstable manifolds: integrate forward
        tspan_fwd = [0 500];
        for s = [-1 1]
            y_initU = Ystar + s*eps_man*v_unst;
            [tU, yU] = ode15s(@(t,y) slow_rhs(t,y,Iapp_pp, ...
                                             gNa,gK,gL,Cm,vNa,vK,vL, ...
                                             m_inf,alpha_n,beta_n), ...
                              tspan_fwd, y_initU);
            if s == -1
                plot(yU(:,1), yU(:,2), 'm', 'LineWidth', 2, ...
                     'DisplayName','unstable manifolds');
            else
                plot(yU(:,1), yU(:,2), 'm', 'LineWidth', 2); % no name
            end
        end

        % Stable manifolds: integrate backward
        tspan_back = [0 -500];
        for s = [-1 1]
            y_initS = Ystar + s*eps_man*v_stab;
            [tS, yS] = ode15s(@(t,y) slow_rhs(t,y,Iapp_pp, ...
                                             gNa,gK,gL,Cm,vNa,vK,vL, ...
                                             m_inf,alpha_n,beta_n), ...
                              tspan_back, y_initS);
            if s == -1
                plot(yS(:,1), yS(:,2), 'c', 'LineWidth', 2, ...
                     'DisplayName','stable manifolds');
            else
                plot(yS(:,1), yS(:,2), 'c', 'LineWidth', 2);
            end
        end
    else
        fprintf('Equilibrium is not a saddle; skipping 1D manifolds.\n');
    end
end

% One legend call, uses DisplayName of all plotted objects
figure(1);
legend('Location','best');


%% 3. Sweep I_app: oscillation threshold and bifurcation diagram
I_vals = linspace(0, 15, 100);
T_long = [0 2000];

Vmax = nan(size(I_vals));
Vmin = nan(size(I_vals));

for k = 1:numel(I_vals)
    Icur = I_vals(k);
    ode_fun_I = @(t,y) slow_rhs(t,y,Icur, ...
                                gNa,gK,gL,Cm,vNa,vK,vL, ...
                                m_inf,alpha_n,beta_n);

    % KEY CHANGE: start from y0, NOT y_cont
    [tL, yL] = ode15s(ode_fun_I, T_long, y0);

    idx    = tL > (T_long(2)/2);
    v_tail = yL(idx,1);

    Vmax(k) = max(v_tail);
    Vmin(k) = min(v_tail);
end

amp    = Vmax - Vmin;
is_osc = amp > 5;  % same threshold if you like


figure(3); clf; hold on;
plot(I_vals(~is_osc), Vmax(~is_osc), 'ko', 'MarkerFaceColor','k'); % steady
plot(I_vals(is_osc),  Vmax(is_osc),  'ro', 'MarkerFaceColor','r'); % upper osc
plot(I_vals(is_osc),  Vmin(is_osc),  'bo', 'MarkerFaceColor','b'); % lower osc
xlabel('I_{app}');
ylabel('v (mV)');
title('Max/min v vs I_{app} (numerical bifurcation diagram)');
legend('equilibrium','oscillation max','oscillation min','Location','best');
box on;

if any(is_osc)
    I_threshold = min(I_vals(is_osc));
    fprintf('Approximate oscillation threshold I_app ≈ %.3f\n', I_threshold);
end

%% Local functions

function [vdot, ndot] = slow_vector_field(V, N, Iapp, ...
    gNa,gK,gL,Cm,vNa,vK,vL,m_inf,alpha_n,beta_n)

    m  = m_inf(V);
    IK = gK .* N.^4 .* (V - vK);
    INa = gNa .* m.^3 .* (0.8 - N) .* (V - vNa);
    IL = gL .* (V - vL);

    vdot = (Iapp - IK - INa - IL) / Cm;

    an = alpha_n(V);
    bn = beta_n(V);
    ndot = an.*(1 - N) - bn.*N;
end

function dy = slow_rhs(~, y, Iapp, ...
    gNa,gK,gL,Cm,vNa,vK,vL,m_inf,alpha_n,beta_n)

    v = y(1);
    n = y(2);

    m   = m_inf(v);
    IK  = gK * n^4           * (v - vK);
    INa = gNa * m^3*(0.8-n)  * (v - vNa);
    IL  = gL * (v - vL);

    dv = (Iapp - IK - INa - IL) / Cm;

    an = alpha_n(v);
    bn = beta_n(v);
    dn = an*(1 - n) - bn*n;

    dy = [dv; dn];
end

function am = alpha_m_fun(v)
    % Safe alpha_m(v) with removable singularity at v = 25 mV
    am = 0.1*(25 - v) ./ (exp((25 - v)/10) - 1);
    idx = abs(25 - v) < 1e-6;
    am(idx) = 1;   % limit at v = 25
end

function an = alpha_n_fun(v)
    % Safe alpha_n(v) with removable singularity at v = 10 mV
    an = 0.01*(10 - v) ./ (exp((10 - v)/10) - 1);
    idx = abs(10 - v) < 1e-6;
    an(idx) = 0.1; % limit at v = 10
end

function J = numerical_jacobian(F, Y)
    % F: R^2 -> R^2, Y: 2x1
    h  = 1e-5;
    F0 = F(Y);
    J  = zeros(2);
    for j = 1:2
        e      = zeros(2,1); e(j) = 1;
        Fp     = F(Y + h*e);
        Fm     = F(Y - h*e);
        J(:,j) = (Fp - Fm)/(2*h);
    end
end

function hh_fast_phase_plane
    % Parameters
    gNa = 120;
    gK  = 36;
    gL  = 0.3;
    Cm  = 1;

    vNa = 115;
    vK  = -12;
    vL  = 10.6;

    h0  = 0.596;
    n0  = 0.3176;

    % Gating-rate functions (v in mV, t in ms)
    alpha_m_raw = @(v) 0.1*(25 - v) ./ (exp((25 - v)/10) - 1);
    beta_m      = @(v) 4*exp(-v/18);

    % Regularize alpha_m at v = 25 mV (remove 0/0)
    alpha_m = @(v) (abs(25-v) < 1e-6).*1 + ...
                   (abs(25-v) >= 1e-6).*alpha_m_raw(v);

    m_inf = @(v) alpha_m(v) ./ (alpha_m(v) + beta_m(v));

    % ODE system: y = [v; m]
    f = @(t,y) [
        ( -gK*n0^4*(y(1)-vK) ...
          -gNa*(y(2)^3)*h0*(y(1)-vNa) ...
          -gL*(y(1)-vL) ) / Cm;
        alpha_m(y(1))*(1 - y(2)) - beta_m(y(1))*y(2)
    ];

    % Phase plane grid
    vmin = -20;  vmax = 120;
    mmin = 0;    mmax = 1;

    [V,M] = meshgrid(linspace(vmin,vmax,25), linspace(mmin,mmax,25));
    dV = (-gK*n0^4.*(V - vK) ...
          -gNa.*(M.^3).*h0.*(V - vNa) ...
          -gL.*(V - vL)) / Cm;
    dM = alpha_m(V).*(1-M) - beta_m(V).*M;

    figure; hold on; box on;
    quiver(V,M,dV,dM,0.4,'Color',[0.7 0.7 0.7]);

    % v-nullcline: dv/dt = 0, v as a function of m
    v_nc_fun = @(m) ( gNa.*(m.^3)*h0*vNa + gK*n0^4*vK + gL*vL ) ./ ...
                    ( gNa.*(m.^3).*h0     + gK*n0^4    + gL );

    m_vals = linspace(mmin,mmax,400);
    v_nc   = v_nc_fun(m_vals);
    plot(v_nc, m_vals, 'r', 'LineWidth', 2);

    % m-nullcline: dm/dt = 0, m = m_inf(v)
    v_vals = linspace(vmin,vmax,400);
    m_nc   = m_inf(v_vals);
    plot(v_vals, m_nc, 'b', 'LineWidth', 2);

    % Find steady states as intersections of nullclines
    F = @(m) m_inf( v_nc_fun(m) ) - m;

    m_scan = linspace(mmin+1e-4,mmax-1e-4,2000);
    F_scan = F(m_scan);
    idx = find(diff(sign(F_scan)));   % sign changes -> roots

    m_eq = zeros(numel(idx),1);
    v_eq = zeros(numel(idx),1);
    for k = 1:numel(idx)
        m_lo = m_scan(idx(k));
        m_hi = m_scan(idx(k)+1);
        m_eq(k) = fzero(F,[m_lo,m_hi]);
        v_eq(k) = v_nc_fun(m_eq(k));
    end

    plot(v_eq, m_eq, 'ko', 'MarkerFaceColor','k', 'MarkerSize',7);

    % Sample trajectories
    tspan = [0 150];  % ms
    ICs = [...
        -10  0.05;
         0   0.2;
        10   0.05;
        20   0.1;
        40   0.15;
        80   0.7;
       100   0.9];

    for k = 1:size(ICs,1)
        [t,y] = ode45(f, tspan, ICs(k,:).');
        plot(y(:,1), y(:,2), 'LineWidth', 1.5);
    end

    xlabel('v (mV)');
    ylabel('m');
    title('Hodgkin-Huxley fast subsystem phase plane');
    legend({'vector field','v-nullcline','m-nullcline','steady states'}, ...
           'Location','bestoutside');
    axis([vmin vmax mmin mmax]);
end

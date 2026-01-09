clear; clc;
g     = 9.81;            % m/s^2
L     = 6.0;             % m
Omega = 2*pi/86400;      % rad/s
w0    = sqrt(g/L);


m     = 28;              
gamma = 0.0;            


h  = 0.05;               % time step (s)
T  = 3*3600;             % 3 hours
t  = (0:h:T)';

% Initial conditions
x0  = 3.0;  y0  = 0.0;
vx0 = 0.0;  vy0 = 0.0;
z0  = [x0; y0; vx0; vy0];

% Angle estimation settings
win_seconds = 600;                      
winN  = max(50, round(win_seconds/h));
stepN = round(winN/4);

warmup_seconds = 900;                   % 15 minutes
warmupN = round(warmup_seconds/stepN/h); % in number of windows (approx)


names   = ["Paris","Durham"];
lat_deg = [48.8566, 54.7761];
alpha_theory_degph = 15*sind(lat_deg);  % deg/hour in 24h model

Omega_hat = zeros(size(lat_deg));
alpha_hat_degph = zeros(size(lat_deg));
phi_hat_deg = zeros(size(lat_deg));

theta_data = cell(1,2);
tc_data    = cell(1,2);


for k = 1:2
    phi   = deg2rad(lat_deg(k));
    sigma = pi/2 - phi;

    Z = rk4_foucault(t, z0, Omega, sigma, w0, gamma, m);
    x = Z(:,1); y = Z(:,2);

    [tc, theta] = estimate_plane_angle_pca(t, x, y, winN, stepN);
    theta_u = unwrap(theta);

    
    idx0 = min(max(1, warmupN), numel(tc)-5);
    tc_fit = tc(idx0:end);
    th_fit = theta_u(idx0:end);

    
    p = polyfit(tc_fit, th_fit, 1);
    Omega_hat(k) = p(1);

    
    alpha_hat_degph(k) = abs(rad2deg(Omega_hat(k))*3600);

    
    ratio = abs(Omega_hat(k))/Omega;
    ratio = max(min(ratio, 1), -1);
    phi_hat_deg(k) = rad2deg(asin(ratio));

    theta_data{k} = theta_u;
    tc_data{k}    = tc;
end


fprintf('\n===== Durham vs Paris results =====\n');
for k = 1:2
    fprintf('%s:\n', names(k));
    fprintf('  true latitude phi = %.4f deg\n', lat_deg(k));
    fprintf('  theory alpha      = %.4f deg/hour\n', alpha_theory_degph(k));
    fprintf('  numeric alpha_hat  = %.4f deg/hour\n', alpha_hat_degph(k));
    fprintf('  inferred phi_hat   = %.4f deg\n\n', phi_hat_deg(k));
end


if ~exist('Plots','dir'); mkdir('Plots'); end


figure; hold on;
for k = 1:2
    plot(tc_data{k}/3600, rad2deg(theta_data{k}), 'LineWidth', 1.4);
end
xlabel('time (hours)','Interpreter','tex');
ylabel('\theta(t) (degrees)','Interpreter','tex');
legend( ...
    sprintf('%s: \\theta(t)', names(1)), ...
    sprintf('%s: \\theta(t)', names(2)), ...
    'Location','best');
grid on;

try
    exportgraphics(gcf, 'Plots/theta_durham_paris.png', 'Resolution', 300);
catch
    print(gcf, 'Plots/theta_durham_paris', '-dpng', '-r300');
end


phi0_deg = lat_deg(1);
phi0 = deg2rad(phi0_deg);
sigma0 = pi/2 - phi0;

T_list = (1:1:10)*3600;
err_deg = zeros(size(T_list));

for i = 1:numel(T_list)
    Ti = T_list(i);
    ti = (0:h:Ti)';

    Z = rk4_foucault(ti, z0, Omega, sigma0, w0, gamma, m);
    x = Z(:,1); y = Z(:,2);

    [tc, theta] = estimate_plane_angle_pca(ti, x, y, winN, stepN);
    theta_u = unwrap(theta);

    % Warm-up trimming for fit
    idx0 = min(max(1, warmupN), numel(tc)-5);
    tc_fit = tc(idx0:end);
    th_fit = theta_u(idx0:end);

    p = polyfit(tc_fit, th_fit, 1);
    Omega_i = p(1);

    ratio = abs(Omega_i)/Omega;
    ratio = max(min(ratio, 1), -1);
    phi_hat = asin(ratio);

    err_deg(i) = abs(rad2deg(phi_hat) - phi0_deg);
end


err_milli = 1e3 * err_deg;

figure;
plot(T_list/3600, err_milli, 'o-', 'LineWidth', 1.2);
xlabel('observation time T (hours)','Interpreter','tex');
ylabel('absolute latitude error $|\hat{\phi}-\phi|\times 10^{-3}$ (degrees)','Interpreter','latex');
grid on;

try
    exportgraphics(gcf, 'Plots/error_vs_time_paris.png', 'Resolution', 300);
catch
    print(gcf, 'Plots/error_vs_time_paris', '-dpng', '-r300');
end
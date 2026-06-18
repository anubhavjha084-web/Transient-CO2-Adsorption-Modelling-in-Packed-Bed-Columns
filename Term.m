% CO₂ Adsorption in a Packed Bed using Finite Differences
clear all; clc; close all;

%% Parameters
L = 1.0;              % Column length (m)
n = 200;             % Number of spatial points
dz = L / (n-1);       % Spatial step

tspan = [0 5000];      % Time span (s)

% System parameters
u = 0.01;             % Superficial velocity (m/s)
D_ax = 1e-6;          % Axial dispersion coefficient (m²/s)
eps = 0.35;            % Bed porosity
rho_s = 50;         % Adsorbent density (kg/m³)
q_max = 1.0;          % Maximum adsorption capacity (mol/kg)
k_a = 0.01;           % Adsorption rate constant (1/s)
k_d = 0.005;           % Desorption rate constant (1/s)

% Initial and boundary conditions
C_in = 1.0;           % Inlet CO₂ concentration (mol/m³)
C0 = zeros(n,1);      % Initial gas concentration
q0 = zeros(n,1);      % Initial adsorbed amount

% Combine into single initial condition vector
y0 = [C0; q0];

%% ODE Solver
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-5);
[t, y] = ode15s(@(t, y) packedBedODE(t, y, n, dz, u, D_ax, eps, rho_s, k_a, k_d, q_max, C_in), tspan, y0, options);

%% Extract solutions
C = y(:,1:n);        % Gas concentration
q = y(:,n+1:end);    % Adsorbed concentration
z = linspace(0, L, n);

%% Plots
figure;

set(gcf); 
set(gca, 'FontSize', 12);

subplot(2,2,1)
plot(z, C(end,:), 'b-', 'LineWidth', 2);
xlabel('Column length z (m)');
ylabel('C(z,t_{end}) [mol/m³]');
title('CO₂ Concentration Profile at Final Time');
grid on;

subplot(2,2,2)
semilogx(t, C(:,end), 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('C(L,t) [mol/m³]');
title('Breakthrough Curve at Column Outlet');
grid on;
ytickformat('%.1e');

subplot(2,2,3)
plot(z, q(end,:), 'g-', 'LineWidth', 2);
xlabel('Column length z (m)');
ylabel('q(z,t_{end}) [mol/kg]');
title('Adsorbed Phase Profile at Final Time');
grid on;

subplot(2,2,4)
semilogx(t, q(:,end), 'm-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('q(L,t) [mol/kg]');
title('Adsorption at Column Outlet Over Time');
grid on;
ytickformat('%.1e');
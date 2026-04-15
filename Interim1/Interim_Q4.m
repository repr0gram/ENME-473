% computes posture, velocity, and acceleration results
% and generates the plots for Question 4
clc;
clear;
close all;

%% Known Values
R1 = sqrt(237.2^2 + 70.6^2);
R2 = 272.4;
R23 = 236.9;
R14 = 279.0;
R3 = 489.9;
R4 = 539.8;
R36 = 342.1;
R46 = 179.3;
R5 = 595.5;
R6 = 378.9;
R7 = 198.2;
R8 = 254.4;
RA = 336.9;         % distance from joint 7/8 to Pin A along link 8
theta1 = deg2rad(180 - atand(70.6/237.2));
alpha = deg2rad(3.2);
beta = deg2rad(3.5);
gamma = deg2rad(1.8);

omega2 = 0.5;   % rad/s
alpha2 = 0;     % rad/s^2

in = 0:120;
theta2_vals = deg2rad(in);

%% Storage
x = deg2rad([140 60 150 150 210 150])';

theta23_vals = zeros(1,length(theta2_vals));
theta14_vals = zeros(1,length(theta2_vals));
theta46_vals = zeros(1,length(theta2_vals));
theta36_vals = zeros(1,length(theta2_vals));
theta8_vals  = zeros(1,length(theta2_vals));
theta7_vals  = zeros(1,length(theta2_vals));

omega3_vals = zeros(1,length(theta2_vals));
omega4_vals = zeros(1,length(theta2_vals));
omega5_vals = zeros(1,length(theta2_vals));
omega6_vals = zeros(1,length(theta2_vals));
omega7_vals = zeros(1,length(theta2_vals));
omega8_vals = zeros(1,length(theta2_vals));

alpha3_vals = zeros(1,length(theta2_vals));
alpha4_vals = zeros(1,length(theta2_vals));
alpha5_vals = zeros(1,length(theta2_vals));
alpha6_vals = zeros(1,length(theta2_vals));
alpha7_vals = zeros(1,length(theta2_vals));
alpha8_vals = zeros(1,length(theta2_vals));

% Q4c: Pin A position, velocity, acceleration storage
Ax = zeros(1,length(theta2_vals));
Ay = zeros(1,length(theta2_vals));
Ax_dot = zeros(1,length(theta2_vals));
Ay_dot = zeros(1,length(theta2_vals));
Ax_ddot = zeros(1,length(theta2_vals));
Ay_ddot = zeros(1,length(theta2_vals));

%% Main Loop: Posture, Velocity, Acceleration
for k = 1:length(theta2_vals)
    theta2 = theta2_vals(k);
    theta23 = x(1);
    theta14 = x(2);
    theta46 = x(3);
    theta36 = x(4);
    theta8  = x(5);
    theta7  = x(6);

    % --- Posture (Newton-Raphson) ---
    while true
        f = [R2*cos(theta2-gamma) + R23*cos(theta23) - R14*cos(theta14) - R1*cos(theta1);
            R2*sin(theta2-gamma) + R23*sin(theta23) - R14*sin(theta14) - R1*sin(theta1);
            R2*cos(theta2-gamma) + R4*cos(theta23) + R46*cos(theta46) - R36*cos(theta36) - R3*cos(theta14+alpha) - R1*cos(theta1);
            R2*sin(theta2-gamma) + R4*sin(theta23) + R46*sin(theta46) - R36*sin(theta36) - R3*sin(theta14+alpha) - R1*sin(theta1);
            R2*cos(theta2-gamma) + R4*cos(theta23) + R6*cos(theta46) - R8*cos(theta8) + R7*cos(theta7) - R5*cos(theta36+beta) - R3*cos(theta14+alpha) - R1*cos(theta1);
            R2*sin(theta2-gamma) + R4*sin(theta23) + R6*sin(theta46) - R8*sin(theta8) + R7*sin(theta7) - R5*sin(theta36+beta) - R3*sin(theta14+alpha) - R1*sin(theta1)];

        J = [-R23*sin(theta23),  R14*sin(theta14),              0,              0,             0,            0;
            R23*cos(theta23), -R14*cos(theta14),              0,              0,             0,            0;
            -R4*sin(theta23),   R3*sin(theta14+alpha), -R46*sin(theta46),  R36*sin(theta36),  0,            0;
            R4*cos(theta23),  -R3*cos(theta14+alpha),  R46*cos(theta46), -R36*cos(theta36),  0,            0;
            -R4*sin(theta23),   R3*sin(theta14+alpha), -R6*sin(theta46),   R5*sin(theta36+beta),  R8*sin(theta8), -R7*sin(theta7);
            R4*cos(theta23),  -R3*cos(theta14+alpha),  R6*cos(theta46),  -R5*cos(theta36+beta), -R8*cos(theta8),  R7*cos(theta7)];

        dx = J\f;
        x_new = x - dx;
        x_new = atan2(sin(x_new), cos(x_new));

        if norm(f,inf) < 1e-6 && norm(x_new - x,inf) < 1e-6
            x = x_new;
            break;
        end
        x = x_new;
        theta23 = x(1); theta14 = x(2); theta46 = x(3);
        theta36 = x(4); theta8 = x(5);  theta7 = x(6);
    end

    % Store posture
    theta23_vals(k) = x(1);
    theta14_vals(k) = x(2);
    theta46_vals(k) = x(3);
    theta36_vals(k) = x(4);
    theta8_vals(k)  = x(5);
    theta7_vals(k)  = x(6);

    % --- Velocity ---
    b_vel = [R2*sin(theta2-gamma)*omega2;
        -R2*cos(theta2-gamma)*omega2;
        R2*sin(theta2-gamma)*omega2;
        -R2*cos(theta2-gamma)*omega2;
        R2*sin(theta2-gamma)*omega2;
        -R2*cos(theta2-gamma)*omega2];

    omega_vec = J \ b_vel;
    % omega_vec = [w23; w14; w46; w36; w8; w7]

    omega3_vals(k) = omega_vec(2);  % w3 = w14
    omega4_vals(k) = omega_vec(1);  % w4 = w23
    omega5_vals(k) = omega_vec(4);  % w5 = w36
    omega6_vals(k) = omega_vec(3);  % w6 = w46
    omega7_vals(k) = omega_vec(6);
    omega8_vals(k) = omega_vec(5);

    % --- Acceleration ---
    w23 = omega_vec(1); w14 = omega_vec(2);
    w46 = omega_vec(3); w36 = omega_vec(4);
    w8  = omega_vec(5); w7  = omega_vec(6);

    rhs_acc = [R2*cos(theta2-gamma)*omega2^2 + R23*cos(theta23)*w23^2 - R14*cos(theta14)*w14^2;
        R2*sin(theta2-gamma)*omega2^2 + R23*sin(theta23)*w23^2 - R14*sin(theta14)*w14^2;
        R2*cos(theta2-gamma)*omega2^2 + R4*cos(theta23)*w23^2 + R46*cos(theta46)*w46^2 - R36*cos(theta36)*w36^2 - R3*cos(theta14+alpha)*w14^2;
        R2*sin(theta2-gamma)*omega2^2 + R4*sin(theta23)*w23^2 + R46*sin(theta46)*w46^2 - R36*sin(theta36)*w36^2 - R3*sin(theta14+alpha)*w14^2;
        R2*cos(theta2-gamma)*omega2^2 + R4*cos(theta23)*w23^2 + R6*cos(theta46)*w46^2 - R8*cos(theta8)*w8^2 + R7*cos(theta7)*w7^2 - R5*cos(theta36+beta)*w36^2 - R3*cos(theta14+alpha)*w14^2;
        R2*sin(theta2-gamma)*omega2^2 + R4*sin(theta23)*w23^2 + R6*sin(theta46)*w46^2 - R8*sin(theta8)*w8^2 + R7*sin(theta7)*w7^2 - R5*sin(theta36+beta)*w36^2 - R3*sin(theta14+alpha)*w14^2];

    alpha_vec = J \ rhs_acc;

    alpha3_vals(k) = alpha_vec(2);
    alpha4_vals(k) = alpha_vec(1);
    alpha5_vals(k) = alpha_vec(4);
    alpha6_vals(k) = alpha_vec(3);
    alpha7_vals(k) = alpha_vec(6);
    alpha8_vals(k) = alpha_vec(5);

    % --- Q4c: Pin A position, velocity, acceleration ---
    a23 = alpha_vec(1);  a46 = alpha_vec(3);  a8 = alpha_vec(5);

    % Position (Origin -> R2 -> R4 -> R6 -> -R8 + RA extension at theta8+deg2rad(168))
    Ax(k) = R2*cos(theta2-gamma) + R4*cos(theta23) + R6*cos(theta46) - R8*cos(theta8) + RA*cos(theta8+deg2rad(168));
    Ay(k) = R2*sin(theta2-gamma) + R4*sin(theta23) + R6*sin(theta46) - R8*sin(theta8) + RA*sin(theta8+deg2rad(168));

    % Velocity (first derivative of position)
    Ax_dot(k) = -R2*sin(theta2-gamma)*omega2 - R4*sin(theta23)*w23 - R6*sin(theta46)*w46 ...
        + R8*sin(theta8)*w8 - RA*sin(theta8+deg2rad(168))*w8;
    Ay_dot(k) =  R2*cos(theta2-gamma)*omega2 + R4*cos(theta23)*w23 + R6*cos(theta46)*w46 ...
        - R8*cos(theta8)*w8 + RA*cos(theta8+deg2rad(168))*w8;

    % Acceleration (second derivative of position, with alpha2 = 0)
    Ax_ddot(k) = -R2*cos(theta2-gamma)*omega2^2 ...
        - R4*cos(theta23)*w23^2 - R4*sin(theta23)*a23 ...
        - R6*cos(theta46)*w46^2 - R6*sin(theta46)*a46 ...
        + R8*cos(theta8)*w8^2 + R8*sin(theta8)*a8 ...
        - RA*cos(theta8+deg2rad(168))*w8^2 - RA*sin(theta8+deg2rad(168))*a8;
    Ay_ddot(k) = -R2*sin(theta2-gamma)*omega2^2 ...
        - R4*sin(theta23)*w23^2 + R4*cos(theta23)*a23 ...
        - R6*sin(theta46)*w46^2 + R6*cos(theta46)*a46 ...
        + R8*sin(theta8)*w8^2 - R8*cos(theta8)*a8 ...
        - RA*sin(theta8+deg2rad(168))*w8^2 + RA*cos(theta8+deg2rad(168))*a8;
end

%% Output directory
% resolve repo-root/exports as the output directory
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = scriptDir;
while ~isfolder(fullfile(repoRoot, '.git'))
    parentDir = fileparts(repoRoot);
    if strcmp(parentDir, repoRoot)
        repoRoot = scriptDir; % fallback if .git not found
        break;
    end
    repoRoot = parentDir;
end
outDir = fullfile(repoRoot, 'exports');
if ~isfolder(outDir); mkdir(outDir); end

%% Plot 1: Angular Velocities
figure
plot(in, omega3_vals, 'LineWidth', 2)
hold on
plot(in, omega4_vals, 'LineWidth', 2)
plot(in, omega5_vals, 'LineWidth', 2)
plot(in, omega6_vals, 'LineWidth', 2)
plot(in, omega7_vals, 'LineWidth', 2)
plot(in, omega8_vals, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Angular Velocity (rad/s)')
title('Angular Velocities vs Input Angle (\omega_2 = 0.5 rad/s)')
legend('\omega_3','\omega_4','\omega_5','\omega_6','\omega_7','\omega_8','Location','best')
grid on
box on
exportgraphics(gcf, fullfile(outDir, 'ENME_473_Project_4a_velocity.png'), 'Resolution', 600);

%% Plot 2: Angular Accelerations
figure
plot(in, alpha3_vals, 'LineWidth', 2)
hold on
plot(in, alpha4_vals, 'LineWidth', 2)
plot(in, alpha5_vals, 'LineWidth', 2)
plot(in, alpha6_vals, 'LineWidth', 2)
plot(in, alpha7_vals, 'LineWidth', 2)
plot(in, alpha8_vals, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Angular Acceleration (rad/s^2)')
title('Angular Accelerations vs Input Angle (\omega_2 = 0.5 rad/s, \alpha_2 = 0)')
legend('\alpha_3','\alpha_4','\alpha_5','\alpha_6','\alpha_7','\alpha_8','Location','best')
grid on
box on
exportgraphics(gcf, fullfile(outDir, 'ENME_473_Project_4a_acceleration.png'), 'Resolution', 600);

%% Write tables to Excel
filename = fullfile(outDir, 'Q4_Results.xlsx');

% Angular Velocities
T_vel = table(in', omega3_vals', omega4_vals', omega5_vals', ...
    omega6_vals', omega7_vals', omega8_vals', ...
    'VariableNames', {'Theta2_deg','omega3','omega4','omega5','omega6','omega7','omega8'});
writetable(T_vel, filename, 'Sheet', 'Angular Velocities');

% Angular Accelerations
T_acc = table(in', alpha3_vals', alpha4_vals', alpha5_vals', ...
    alpha6_vals', alpha7_vals', alpha8_vals', ...
    'VariableNames', {'Theta2_deg','alpha3','alpha4','alpha5','alpha6','alpha7','alpha8'});
writetable(T_acc, filename, 'Sheet', 'Angular Accelerations');

fprintf('Q4 results written to %s\n', filename);

%% Print summary tables for key angles
idx = [1 31 61 91 121]; % theta2 = 0, 30, 60, 90, 120

fprintf('\n=== Angular Velocities (rad/s) ===\n');
fprintf('%10s %10s %10s %10s %10s %10s %10s\n', ...
    'Theta 2', 'omega3', 'omega4', 'omega5', 'omega6', 'omega7', 'omega8');
for k = idx
    fprintf('%10d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', ...
        in(k), omega3_vals(k), omega4_vals(k), omega5_vals(k), ...
        omega6_vals(k), omega7_vals(k), omega8_vals(k));
end

fprintf('\n=== Angular Accelerations (rad/s^2) ===\n');
fprintf('%10s %10s %10s %10s %10s %10s %10s\n', ...
    'Theta 2', 'alpha3', 'alpha4', 'alpha5', 'alpha6', 'alpha7', 'alpha8');
for k = idx
    fprintf('%10d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', ...
        in(k), alpha3_vals(k), alpha4_vals(k), alpha5_vals(k), ...
        alpha6_vals(k), alpha7_vals(k), alpha8_vals(k));
end

%% Q4c: Pin A Plots

% Plot 3: Pin A Horizontal and Vertical Acceleration
figure
plot(in, Ax_ddot, 'LineWidth', 2)
hold on
plot(in, Ay_ddot, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Acceleration (mm/s^2)')
title('Pin A Linear Acceleration vs Input Angle (\omega_2 = 0.5 rad/s, \alpha_2 = 0)')
legend('a_{Ax} (horizontal)', 'a_{Ay} (vertical)', 'Location', 'best')
grid on
box on
exportgraphics(gcf, fullfile(outDir, 'ENME_473_Project_4c_acceleration.png'), 'Resolution', 600);

% Plot 4: Pin A Horizontal and Vertical Velocity
figure
plot(in, Ax_dot, 'LineWidth', 2)
hold on
plot(in, Ay_dot, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Velocity (mm/s)')
title('Pin A Linear Velocity vs Input Angle (\omega_2 = 0.5 rad/s)')
legend('v_{Ax} (horizontal)', 'v_{Ay} (vertical)', 'Location', 'best')
grid on
box on
exportgraphics(gcf, fullfile(outDir, 'ENME_473_Project_4c_velocity.png'), 'Resolution', 600);

%% Q4c: Pin A Tables

% Write to Excel
T_pinA = table(in', Ax', Ay', Ax_dot', Ay_dot', Ax_ddot', Ay_ddot', ...
    'VariableNames', {'Theta2_deg','X_A_mm','Y_A_mm','Vx_A_mm_s','Vy_A_mm_s','Ax_A_mm_s2','Ay_A_mm_s2'});
writetable(T_pinA, filename, 'Sheet', 'Pin A');

% Print summary
fprintf('\n=== Pin A Velocity (mm/s) ===\n');
fprintf('%10s %12s %12s\n', 'Theta 2', 'Vx_A', 'Vy_A');
for k = idx
    fprintf('%10d %12.4f %12.4f\n', in(k), Ax_dot(k), Ay_dot(k));
end

fprintf('\n=== Pin A Acceleration (mm/s^2) ===\n');
fprintf('%10s %12s %12s\n', 'Theta 2', 'Ax_A', 'Ay_A');
for k = idx
    fprintf('%10d %12.4f %12.4f\n', in(k), Ax_ddot(k), Ay_ddot(k));
end

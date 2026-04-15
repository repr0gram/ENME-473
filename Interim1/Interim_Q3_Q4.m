% ENME 473 Project Deliverable 1 - Questions 3 and 4
% Posture, velocity, acceleration, and Pin A kinematics for the
% convertible roof mechanism over theta2 = 0 to 120 deg.
clc; clear; close all;

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
RB = 580.7;         % third side of the R8-RA-RB triangle at joint 7/8
theta1 = deg2rad(180 - atand(70.6/237.2));
alpha = deg2rad(3.2);
beta = deg2rad(3.5);
gamma = deg2rad(1.8);

% angle between R8 and RA at the joint that hosts pin A (law of cosines)
PhiA = acos((RA^2 + R8^2 - RB^2) / (2*RA*R8));

omega2 = 0.1745;   % rad/s (input angular velocity)
alpha2 = 0;     % rad/s^2 (input angular acceleration)

in = 0:120;
theta2_vals = deg2rad(in);
N = length(theta2_vals);

%% Storage
% initial guesses for unknown link angles (read off reference at theta2 = 0)
% order: [theta23, theta14, theta46, theta36, theta8, theta7]
% theta23 = theta4 - gamma ~ 168, theta14 = theta3 - alpha ~ -3,
% theta46 = theta6 ~ 10, theta36 = theta5 - beta ~ 167,
% theta8_solved = theta8_displayed - 180 ~ -2, theta7_solved = theta7_displayed - 180 ~ -170
x = deg2rad([168 -3 10 167 -2 -170])';

% posture
theta23_vals = zeros(1, N);
theta14_vals = zeros(1, N);
theta46_vals = zeros(1, N);
theta36_vals = zeros(1, N);
theta8_vals  = zeros(1, N);
theta7_vals  = zeros(1, N);

% angular velocities
omega3_vals = zeros(1, N);
omega4_vals = zeros(1, N);
omega5_vals = zeros(1, N);
omega6_vals = zeros(1, N);
omega7_vals = zeros(1, N);
omega8_vals = zeros(1, N);

% angular accelerations
alpha3_vals = zeros(1, N);
alpha4_vals = zeros(1, N);
alpha5_vals = zeros(1, N);
alpha6_vals = zeros(1, N);
alpha7_vals = zeros(1, N);
alpha8_vals = zeros(1, N);

% Pin A position, velocity, acceleration
Ax = zeros(1, N);
Ay = zeros(1, N);
Ax_dot = zeros(1, N);
Ay_dot = zeros(1, N);
Ax_ddot = zeros(1, N);
Ay_ddot = zeros(1, N);

%% Main Loop: posture, velocity, acceleration, Pin A kinematics
for k = 1:N
    theta2 = theta2_vals(k);
    theta23 = x(1);
    theta14 = x(2);
    theta46 = x(3);
    theta36 = x(4);
    theta8  = x(5);
    theta7  = x(6);

    % --- Posture (Newton-Raphson) ---
    maxIter = 50;
    iter = 0;
    while true
        iter = iter + 1;
        f = [R2*cos(theta2) + R23*cos(theta23) - R14*cos(theta14) - R1*cos(theta1);
            R2*sin(theta2) + R23*sin(theta23) - R14*sin(theta14) - R1*sin(theta1);
            R2*cos(theta2) + R4*cos(theta23+gamma) + R46*cos(theta46) - R36*cos(theta36) - R3*cos(theta14+alpha) - R1*cos(theta1);
            R2*sin(theta2) + R4*sin(theta23+gamma) + R46*sin(theta46) - R36*sin(theta36) - R3*sin(theta14+alpha) - R1*sin(theta1);
            R2*cos(theta2) + R4*cos(theta23+gamma) + R6*cos(theta46) - R8*cos(theta8) + R7*cos(theta7) - R5*cos(theta36+beta) - R3*cos(theta14+alpha) - R1*cos(theta1);
            R2*sin(theta2) + R4*sin(theta23+gamma) + R6*sin(theta46) - R8*sin(theta8) + R7*sin(theta7) - R5*sin(theta36+beta) - R3*sin(theta14+alpha) - R1*sin(theta1)];

        J = [-R23*sin(theta23),        R14*sin(theta14),              0,              0,             0,            0;
              R23*cos(theta23),       -R14*cos(theta14),              0,              0,             0,            0;
             -R4*sin(theta23+gamma),   R3*sin(theta14+alpha), -R46*sin(theta46),  R36*sin(theta36),  0,            0;
              R4*cos(theta23+gamma),  -R3*cos(theta14+alpha),  R46*cos(theta46), -R36*cos(theta36),  0,            0;
             -R4*sin(theta23+gamma),   R3*sin(theta14+alpha), -R6*sin(theta46),   R5*sin(theta36+beta),  R8*sin(theta8), -R7*sin(theta7);
              R4*cos(theta23+gamma),  -R3*cos(theta14+alpha),  R6*cos(theta46),  -R5*cos(theta36+beta), -R8*cos(theta8),  R7*cos(theta7)];

        dx = J\f;
        x_new = x - dx;
        x_new = atan2(sin(x_new), cos(x_new));

        if norm(f,inf) < 1e-6 && norm(x_new - x,inf) < 1e-6
            x = x_new;
            break;
        end
        if iter >= maxIter
            warning('Newton-Raphson did not converge at theta2 = %d deg', in(k));
            x = x_new;
            break;
        end
        x = x_new;
        theta23 = x(1); theta14 = x(2); theta46 = x(3);
        theta36 = x(4); theta8 = x(5);  theta7 = x(6);
    end

    theta23_vals(k) = x(1);
    theta14_vals(k) = x(2);
    theta46_vals(k) = x(3);
    theta36_vals(k) = x(4);
    theta8_vals(k)  = x(5);
    theta7_vals(k)  = x(6);

    % --- Jacobian conditioning check ---
    % warn if the converged J is near-singular (velocity/accel will be unreliable)
    rc = rcond(J);
    if rc < 1e-12
        warning('Near-singular Jacobian at theta2 = %d deg: rcond = %.3e', in(k), rc);
    end

    % --- Velocity ---
    b_vel = [R2*sin(theta2)*omega2;
            -R2*cos(theta2)*omega2;
             R2*sin(theta2)*omega2;
            -R2*cos(theta2)*omega2;
             R2*sin(theta2)*omega2;
            -R2*cos(theta2)*omega2];

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

    rhs_acc = [R2*cos(theta2)*omega2^2 + R23*cos(theta23)*w23^2 - R14*cos(theta14)*w14^2;
        R2*sin(theta2)*omega2^2 + R23*sin(theta23)*w23^2 - R14*sin(theta14)*w14^2;
        R2*cos(theta2)*omega2^2 + R4*cos(theta23+gamma)*w23^2 + R46*cos(theta46)*w46^2 - R36*cos(theta36)*w36^2 - R3*cos(theta14+alpha)*w14^2;
        R2*sin(theta2)*omega2^2 + R4*sin(theta23+gamma)*w23^2 + R46*sin(theta46)*w46^2 - R36*sin(theta36)*w36^2 - R3*sin(theta14+alpha)*w14^2;
        R2*cos(theta2)*omega2^2 + R4*cos(theta23+gamma)*w23^2 + R6*cos(theta46)*w46^2 - R8*cos(theta8)*w8^2 + R7*cos(theta7)*w7^2 - R5*cos(theta36+beta)*w36^2 - R3*cos(theta14+alpha)*w14^2;
        R2*sin(theta2)*omega2^2 + R4*sin(theta23+gamma)*w23^2 + R6*sin(theta46)*w46^2 - R8*sin(theta8)*w8^2 + R7*sin(theta7)*w7^2 - R5*sin(theta36+beta)*w36^2 - R3*sin(theta14+alpha)*w14^2];

    alpha_vec = J \ rhs_acc;

    alpha3_vals(k) = alpha_vec(2);
    alpha4_vals(k) = alpha_vec(1);
    alpha5_vals(k) = alpha_vec(4);
    alpha6_vals(k) = alpha_vec(3);
    alpha7_vals(k) = alpha_vec(6);
    alpha8_vals(k) = alpha_vec(5);

    % --- Pin A kinematics ---
    a23 = alpha_vec(1);  a46 = alpha_vec(3);  a8 = alpha_vec(5);

    Ax(k) = R2*cos(theta2) + R4*cos(theta23+gamma) + R6*cos(theta46) - R8*cos(theta8) + RA*cos(theta8+PhiA);
    Ay(k) = R2*sin(theta2) + R4*sin(theta23+gamma) + R6*sin(theta46) - R8*sin(theta8) + RA*sin(theta8+PhiA);

    Ax_dot(k) = -R2*sin(theta2)*omega2 - R4*sin(theta23+gamma)*w23 - R6*sin(theta46)*w46 ...
        + R8*sin(theta8)*w8 - RA*sin(theta8+PhiA)*w8;
    Ay_dot(k) =  R2*cos(theta2)*omega2 + R4*cos(theta23+gamma)*w23 + R6*cos(theta46)*w46 ...
        - R8*cos(theta8)*w8 + RA*cos(theta8+PhiA)*w8;

    Ax_ddot(k) = -R2*cos(theta2)*omega2^2 ...
        - R4*cos(theta23+gamma)*w23^2 - R4*sin(theta23+gamma)*a23 ...
        - R6*cos(theta46)*w46^2 - R6*sin(theta46)*a46 ...
        + R8*cos(theta8)*w8^2 + R8*sin(theta8)*a8 ...
        - RA*cos(theta8+PhiA)*w8^2 - RA*sin(theta8+PhiA)*a8;
    Ay_ddot(k) = -R2*sin(theta2)*omega2^2 ...
        - R4*sin(theta23+gamma)*w23^2 + R4*cos(theta23+gamma)*a23 ...
        - R6*sin(theta46)*w46^2 + R6*cos(theta46)*a46 ...
        + R8*sin(theta8)*w8^2 - R8*cos(theta8)*a8 ...
        - RA*sin(theta8+PhiA)*w8^2 + RA*cos(theta8+PhiA)*a8;
end

%% Link angles (degrees) for Q3a plot/table
theta23_deg = rad2deg(unwrap(theta23_vals));
theta14_deg = rad2deg(unwrap(theta14_vals));
theta46_deg = rad2deg(unwrap(theta46_vals));
theta36_deg = rad2deg(unwrap(theta36_vals));
theta8_deg  = rad2deg(unwrap(theta8_vals));
theta7_deg  = rad2deg(unwrap(theta7_vals));

theta3_deg = theta14_deg + rad2deg(alpha);
theta4_deg = theta23_deg + rad2deg(gamma);
theta5_deg = theta36_deg + rad2deg(beta);
theta6_deg = theta46_deg;
theta7_deg = theta7_deg + 180;
theta8_deg = theta8_deg + 180;

%% Kinematic coefficients (Q4a)
% first-order:  h_i  = d(theta_i)/d(theta_2) = omega_i / omega_2
% second-order: h_i' = d^2(theta_i)/d(theta_2)^2 = (alpha_i - h_i*alpha_2) / omega_2^2
h3  = omega3_vals / omega2;
h4  = omega4_vals / omega2;
h5  = omega5_vals / omega2;
h6  = omega6_vals / omega2;
h7  = omega7_vals / omega2;
h8  = omega8_vals / omega2;

hp3 = (alpha3_vals - h3*alpha2) / omega2^2;
hp4 = (alpha4_vals - h4*alpha2) / omega2^2;
hp5 = (alpha5_vals - h5*alpha2) / omega2^2;
hp6 = (alpha6_vals - h6*alpha2) / omega2^2;
hp7 = (alpha7_vals - h7*alpha2) / omega2^2;
hp8 = (alpha8_vals - h8*alpha2) / omega2^2;

% Pin A point-path coefficients (mm/rad and mm/rad^2)
Ax_h  = Ax_dot  / omega2;
Ay_h  = Ay_dot  / omega2;
Ax_hp = (Ax_ddot - Ax_h*alpha2) / omega2^2;
Ay_hp = (Ay_ddot - Ay_h*alpha2) / omega2^2;

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

%% Plot (Q3a): Link Angles
figure
plot(in, theta3_deg, 'LineWidth', 2); hold on
plot(in, theta4_deg, 'LineWidth', 2)
plot(in, theta5_deg, 'LineWidth', 2)
plot(in, theta6_deg, 'LineWidth', 2)
plot(in, theta7_deg, 'LineWidth', 2)
plot(in, theta8_deg, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Link Angle (deg)')
title('Link Angles vs Input Angle')
legend('\theta_3','\theta_4','\theta_5','\theta_6','\theta_7','\theta_8','Location','best')
grid on; box on
exportgraphics(gcf, fullfile(outDir, 'ENME_473_Project_3a.png'), 'Resolution', 600);

%% Plot (Q3c): Pin A Path
figure
plot(Ax, Ay, 'b-', 'LineWidth', 2); hold on
plot(Ax(1),   Ay(1),   'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', '\theta_2 = 0°')
plot(Ax(end), Ay(end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', '\theta_2 = 120°')
for idx = [31 61 91]
    plot(Ax(idx), Ay(idx), 'kd', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off')
    text(Ax(idx)+10, Ay(idx)+10, sprintf('\\theta_2 = %d°', in(idx)), 'FontSize', 9)
end
xlabel('X (mm)')
ylabel('Y (mm)')
title('Path Traced by Pin A')
legend('Path of Pin A', '\theta_2 = 0°', '\theta_2 = 120°', 'Location', 'best')
grid on; box on; axis equal
exportgraphics(gcf, fullfile(outDir, 'ENME_473_Project_3c.png'), 'Resolution', 600);

%% Plot (Q4a): First-order Kinematic Coefficients
figure
plot(in, h3, 'LineWidth', 2); hold on
plot(in, h4, 'LineWidth', 2)
plot(in, h5, 'LineWidth', 2)
plot(in, h6, 'LineWidth', 2)
plot(in, h7, 'LineWidth', 2)
plot(in, h8, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('First-order kinematic coefficient')
title('First-order Kinematic Coefficients vs Input Angle')
legend('kinC3','kinC4','kinC5','kinC6','kinC7','kinC8','Location','best')
grid on; box on
exportgraphics(gcf, fullfile(outDir, 'ENME_473_Project_4a_first_order.png'), 'Resolution', 600);

%% Plot (Q4a): Second-order Kinematic Coefficients
figure
plot(in, hp3, 'LineWidth', 2); hold on
plot(in, hp4, 'LineWidth', 2)
plot(in, hp5, 'LineWidth', 2)
plot(in, hp6, 'LineWidth', 2)
plot(in, hp7, 'LineWidth', 2)
plot(in, hp8, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Second-order kinematic coefficient')
title('Second-order Kinematic Coefficients vs Input Angle')
legend('kinC3''','kinC4''','kinC5''','kinC6''','kinC7''','kinC8''','Location','best')
grid on; box on
exportgraphics(gcf, fullfile(outDir, 'ENME_473_Project_4a_second_order.png'), 'Resolution', 600);

%% Plot (Q4a): Angular Velocities
figure
plot(in, rad2deg(omega3_vals), 'LineWidth', 2); hold on
plot(in, rad2deg(omega4_vals), 'LineWidth', 2)
plot(in, rad2deg(omega5_vals), 'LineWidth', 2)
plot(in, rad2deg(omega6_vals), 'LineWidth', 2)
plot(in, rad2deg(omega7_vals), 'LineWidth', 2)
plot(in, rad2deg(omega8_vals), 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Angular velocity (deg/s)')
title('Angular Velocities of All Links vs Input Angle')
legend('\omega_3','\omega_4','\omega_5','\omega_6','\omega_7','\omega_8','Location','best')
grid on; box on
exportgraphics(gcf, fullfile(outDir, 'ENME_473_Project_4a_velocity.png'), 'Resolution', 600);

%% Plot (Q4b): Angular Accelerations
figure
plot(in, rad2deg(alpha3_vals), 'LineWidth', 2); hold on
plot(in, rad2deg(alpha4_vals), 'LineWidth', 2)
plot(in, rad2deg(alpha5_vals), 'LineWidth', 2)
plot(in, rad2deg(alpha6_vals), 'LineWidth', 2)
plot(in, rad2deg(alpha7_vals), 'LineWidth', 2)
plot(in, rad2deg(alpha8_vals), 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Angular acceleration (deg/s^2)')
title('Angular Accelerations of All Links vs Input Angle')
legend('\alpha_3','\alpha_4','\alpha_5','\alpha_6','\alpha_7','\alpha_8','Location','best')
grid on; box on
exportgraphics(gcf, fullfile(outDir, 'ENME_473_Project_4a_acceleration.png'), 'Resolution', 600);

%% Plot (Q4c): Pin A Linear Velocity
figure
plot(in, Ax_dot, 'LineWidth', 2); hold on
plot(in, Ay_dot, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Velocity of Pin A (mm/s)')
title('XY Velocity of Pin A vs Input Angle')
legend('V_x', 'V_y', 'Location', 'best')
grid on; box on
exportgraphics(gcf, fullfile(outDir, 'ENME_473_Project_4c_velocity.png'), 'Resolution', 600);

%% Plot (Q4c): Pin A Linear Acceleration
figure
plot(in, Ax_ddot, 'LineWidth', 2); hold on
plot(in, Ay_ddot, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Acceleration of Pin A (mm/s^2)')
title('XY Acceleration of Pin A vs Input Angle')
legend('A_x', 'A_y', 'Location', 'best')
grid on; box on
exportgraphics(gcf, fullfile(outDir, 'ENME_473_Project_4c_acceleration.png'), 'Resolution', 600);

%% Excel output
filename = fullfile(outDir, 'Q3_Q4_Results.xlsx');

T_angles = table(in', theta3_deg', theta4_deg', theta5_deg', ...
    theta6_deg', theta7_deg', theta8_deg', ...
    'VariableNames', {'Theta2_deg','Theta3','Theta4','Theta5','Theta6','Theta7','Theta8'});
writetable(T_angles, filename, 'Sheet', 'Link Angles');

T_h = table(in', h3', h4', h5', h6', h7', h8', ...
    'VariableNames', {'Theta2_deg','h3','h4','h5','h6','h7','h8'});
writetable(T_h, filename, 'Sheet', 'First-order Coefficients');

T_hp = table(in', hp3', hp4', hp5', hp6', hp7', hp8', ...
    'VariableNames', {'Theta2_deg','hp3','hp4','hp5','hp6','hp7','hp8'});
writetable(T_hp, filename, 'Sheet', 'Second-order Coefficients');

T_vel = table(in', omega3_vals', omega4_vals', omega5_vals', ...
    omega6_vals', omega7_vals', omega8_vals', ...
    'VariableNames', {'Theta2_deg','omega3','omega4','omega5','omega6','omega7','omega8'});
writetable(T_vel, filename, 'Sheet', 'Angular Velocities');

T_acc = table(in', alpha3_vals', alpha4_vals', alpha5_vals', ...
    alpha6_vals', alpha7_vals', alpha8_vals', ...
    'VariableNames', {'Theta2_deg','alpha3','alpha4','alpha5','alpha6','alpha7','alpha8'});
writetable(T_acc, filename, 'Sheet', 'Angular Accelerations');

T_pinA = table(in', Ax', Ay', Ax_dot', Ay_dot', Ax_ddot', Ay_ddot', ...
    'VariableNames', {'Theta2_deg','X_A_mm','Y_A_mm','Vx_A_mm_s','Vy_A_mm_s','Ax_A_mm_s2','Ay_A_mm_s2'});
writetable(T_pinA, filename, 'Sheet', 'Pin A Kinematics');

fprintf('Q3 & Q4 results written to %s\n', filename);

%% Summary tables for key input angles
idx = [1 31 61 91 121]; % theta2 = 0, 30, 60, 90, 120

fprintf('\n=== Link Angles (degrees) ===\n');
fprintf('%10s %10s %10s %10s %10s %10s %10s\n', ...
    'Theta 2', 'Theta 3', 'Theta 4', 'Theta 5', 'Theta 6', 'Theta 7', 'Theta 8');
for k = idx
    fprintf('%10d %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n', ...
        in(k), theta3_deg(k), theta4_deg(k), theta5_deg(k), ...
        theta6_deg(k), theta7_deg(k), theta8_deg(k));
end

fprintf('\n=== First-order Kinematic Coefficients (h_i) ===\n');
fprintf('%10s %10s %10s %10s %10s %10s %10s\n', ...
    'Theta 2', 'h3', 'h4', 'h5', 'h6', 'h7', 'h8');
for k = idx
    fprintf('%10d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', ...
        in(k), h3(k), h4(k), h5(k), h6(k), h7(k), h8(k));
end

fprintf('\n=== Second-order Kinematic Coefficients (h_i'') ===\n');
fprintf('%10s %10s %10s %10s %10s %10s %10s\n', ...
    'Theta 2', 'hp3', 'hp4', 'hp5', 'hp6', 'hp7', 'hp8');
for k = idx
    fprintf('%10d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', ...
        in(k), hp3(k), hp4(k), hp5(k), hp6(k), hp7(k), hp8(k));
end

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

fprintf('\n=== Pin A Position (mm) ===\n');
fprintf('%10s %12s %12s\n', 'Theta 2', 'X_A', 'Y_A');
for k = idx
    fprintf('%10d %12.2f %12.2f\n', in(k), Ax(k), Ay(k));
end

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

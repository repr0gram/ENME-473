% ENME 473 Project Deliverable 1 - Questions 3a & 3c
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
theta1 = deg2rad(180 - atand(70.6/237.2));
alpha = deg2rad(3.2);
beta = deg2rad(3.5);
gamma = deg2rad(1.8);

% create the range of input angles
in = 0:120;
theta2_vals = deg2rad(in);

%% Posture Analysis (Newton-Raphson)

% store initial guesses
x = deg2rad([140 60 150 150 210 150])';

% create storage for angle values (radians)
theta23_vals = zeros(1, length(theta2_vals));
theta14_vals = zeros(1, length(theta2_vals));
theta46_vals = zeros(1, length(theta2_vals));
theta36_vals = zeros(1, length(theta2_vals));
theta8_vals  = zeros(1, length(theta2_vals));
theta7_vals  = zeros(1, length(theta2_vals));

% storage for Pin A coordinates
Ax = zeros(1, length(theta2_vals));
Ay = zeros(1, length(theta2_vals));

% storage for angular velocities
omegas = zeros(6, length(theta2_vals));
omega2 = 0.5; % input angular velocity (rad/s)

% loop through input values
for k = 1:length(theta2_vals)
    theta2 = theta2_vals(k);
    theta23 = x(1); theta14 = x(2); theta46 = x(3);
    theta36 = x(4); theta8  = x(5); theta7  = x(6);

    while true
        % store the force matrix
        f = [R2*cos(theta2-gamma) + R23*cos(theta23) - R14*cos(theta14) - R1*cos(theta1);
            R2*sin(theta2-gamma) + R23*sin(theta23) - R14*sin(theta14) - R1*sin(theta1);
            R2*cos(theta2-gamma) + R4*cos(theta23) + R46*cos(theta46) - R36*cos(theta36) - R3*cos(theta14+alpha) - R1*cos(theta1);
            R2*sin(theta2-gamma) + R4*sin(theta23) + R46*sin(theta46) - R36*sin(theta36) - R3*sin(theta14+alpha) - R1*sin(theta1);
            R2*cos(theta2-gamma) + R4*cos(theta23) + R6*cos(theta46) - R8*cos(theta8) + R7*cos(theta7) - R5*cos(theta36+beta) - R3*cos(theta14+alpha) - R1*cos(theta1);
            R2*sin(theta2-gamma) + R4*sin(theta23) + R6*sin(theta46) - R8*sin(theta8) + R7*sin(theta7) - R5*sin(theta36+beta) - R3*sin(theta14+alpha) - R1*sin(theta1)];

        % store the jacobian
        J = [-R23*sin(theta23),  R14*sin(theta14),              0,                    0,             0,            0;
            R23*cos(theta23), -R14*cos(theta14),              0,                    0,             0,            0;
            -R4*sin(theta23),   R3*sin(theta14+alpha), -R46*sin(theta46),  R36*sin(theta36),       0,            0;
            R4*cos(theta23),  -R3*cos(theta14+alpha),  R46*cos(theta46), -R36*cos(theta36),       0,            0;
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
        theta36 = x(4); theta8  = x(5); theta7  = x(6);
    end

    % store posture results
    theta23_vals(k) = x(1);
    theta14_vals(k) = x(2);
    theta46_vals(k) = x(3);
    theta36_vals(k) = x(4);
    theta8_vals(k)  = x(5);
    theta7_vals(k)  = x(6);

    % velocity analysis: J * omega = -df/dtheta2 * omega2
    theta23 = x(1); theta14 = x(2); theta46 = x(3);
    theta36 = x(4); theta8  = x(5); theta7  = x(6);
    J = [-R23*sin(theta23),  R14*sin(theta14),              0,                    0,             0,            0;
          R23*cos(theta23), -R14*cos(theta14),              0,                    0,             0,            0;
         -R4*sin(theta23),   R3*sin(theta14+alpha), -R46*sin(theta46),  R36*sin(theta36),       0,            0;
          R4*cos(theta23),  -R3*cos(theta14+alpha),  R46*cos(theta46), -R36*cos(theta36),       0,            0;
         -R4*sin(theta23),   R3*sin(theta14+alpha), -R6*sin(theta46),   R5*sin(theta36+beta),  R8*sin(theta8), -R7*sin(theta7);
          R4*cos(theta23),  -R3*cos(theta14+alpha),  R6*cos(theta46),  -R5*cos(theta36+beta), -R8*cos(theta8),  R7*cos(theta7)];
    df_dt2 = [-R2*sin(theta2-gamma);
               R2*cos(theta2-gamma);
              -R2*sin(theta2-gamma);
               R2*cos(theta2-gamma);
              -R2*sin(theta2-gamma);
               R2*cos(theta2-gamma)];
    omegas(:,k) = J \ (-df_dt2 * omega2);

    % compute Pin A position: Origin -> R2 -> R4 -> R6 -> (R8+RA) extension
    % theta4 = theta23, theta6 = theta46
    Ax(k) = R2*cos(theta2-gamma) + R4*cos(x(1)) + R6*cos(x(3)) - R8*cos(x(5)) + RA*cos(x(5)+168);
    Ay(k) = R2*sin(theta2-gamma) + R4*sin(x(1)) + R6*sin(x(3)) - R8*sin(x(5)) + RA*sin(x(5)+168);
end

%% Post-process Angles

% convert to degrees with unwrapping
theta23_vals = rad2deg(unwrap(theta23_vals));
theta14_vals = rad2deg(unwrap(theta14_vals));
theta46_vals = rad2deg(unwrap(theta46_vals));
theta36_vals = rad2deg(unwrap(theta36_vals));
theta8_vals  = rad2deg(unwrap(theta8_vals));
theta7_vals  = rad2deg(unwrap(theta7_vals));

% calculate actual link angles using constraint equations
theta3_vals = theta14_vals + rad2deg(alpha);
theta4_vals = theta23_vals + rad2deg(gamma);
theta5_vals = theta36_vals + rad2deg(beta);
theta6_vals = theta46_vals;
theta7_vals = theta7_vals +180;
theta8_vals = theta8_vals +180;

outDir = fileparts(mfilename('fullpath'));

%% Plot 1: Link Angles vs Input Angle (Q3a)
figure
plot(in, theta3_vals, 'LineWidth', 2)
hold on
plot(in, theta4_vals, 'LineWidth', 2)
plot(in, theta5_vals, 'LineWidth', 2)
plot(in, theta6_vals, 'LineWidth', 2)
plot(in, theta7_vals, 'LineWidth', 2)
plot(in, theta8_vals, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Link Angle (deg)')
title('Link Angles vs Input Angle')
legend('\theta_3','\theta_4','\theta_5','\theta_6','\theta_7','\theta_8','Location','best')
grid on
box on
exportgraphics(gcf, fullfile(outDir, 'ENME_473_Project_3a.png'), 'Resolution', 600);

%% Plot 2: Path Traced by Pin A (Q3c)
figure
plot(Ax, Ay, 'b-', 'LineWidth', 2)
hold on
% mark start and end points
plot(Ax(1),   Ay(1),   'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', '\theta_2 = 0°')
plot(Ax(end), Ay(end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', '\theta_2 = 120°')
% mark every 30 degrees
for idx = [31 61 91]  % theta2 = 30, 60, 90
    plot(Ax(idx), Ay(idx), 'kd', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off')
    text(Ax(idx)+10, Ay(idx)+10, sprintf('\\theta_2 = %d°', in(idx)), 'FontSize', 9)
end
xlabel('X (mm)')
ylabel('Y (mm)')
title('Path Traced by Pin A')
legend('Path of Pin A', '\theta_2 = 0°', '\theta_2 = 120°', 'Location', 'best')
grid on
box on
axis equal
exportgraphics(gcf, fullfile(outDir, 'ENME_473_Project_3c.png'), 'Resolution', 600);

%% Write tables to Excel
filename = fullfile(outDir, 'Q3_Results.xlsx');

% Q3a - Link Angles
T_angles = table(in', theta3_vals', theta4_vals', theta5_vals', ...
    theta6_vals', theta7_vals', theta8_vals', ...
    'VariableNames', {'Theta2_deg','Theta3','Theta4','Theta5','Theta6','Theta7','Theta8'});
writetable(T_angles, filename, 'Sheet', 'Link Angles');

% Q3b - Angular Velocities
% omegas rows: [t23, t14, t46, t36, t8, t7]
T_vel = table(in', omegas(2,:)', omegas(1,:)', omegas(4,:)', ...
    omegas(3,:)', omegas(6,:)', omegas(5,:)', ...
    'VariableNames', {'Theta2_deg','omega3','omega4','omega5','omega6','omega7','omega8'});
writetable(T_vel, filename, 'Sheet', 'Angular Velocities');

% Q3c - Pin A Position
T_pinA = table(in', Ax', Ay', ...
    'VariableNames', {'Theta2_deg','X_A_mm','Y_A_mm'});
writetable(T_pinA, filename, 'Sheet', 'Pin A Position');

fprintf('Q3 results written to %s\n', filename);

%% Print summary tables for key angles
idx = [1 31 61 91 121]; % theta2 = 0, 30, 60, 90, 120

fprintf('\n=== Link Angles (degrees) ===\n');
fprintf('%10s %10s %10s %10s %10s %10s %10s\n', ...
    'Theta 2', 'Theta 3', 'Theta 4', 'Theta 5', 'Theta 6', 'Theta 7', 'Theta 8');
for k = idx
    fprintf('%10d %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n', ...
        in(k), theta3_vals(k), theta4_vals(k), theta5_vals(k), ...
        theta6_vals(k), theta7_vals(k), theta8_vals(k));
end

fprintf('\n=== Angular Velocities (rad/s) ===\n');
fprintf('%10s %10s %10s %10s %10s %10s %10s\n', ...
    'Theta 2', 'omega3', 'omega4', 'omega5', 'omega6', 'omega7', 'omega8');
for k = idx
    fprintf('%10d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', ...
        in(k), omegas(2,k), omegas(1,k), omegas(4,k), ...
        omegas(3,k), omegas(6,k), omegas(5,k));
end

fprintf('\n=== Pin A Position (mm) ===\n');
fprintf('%10s %12s %12s\n', 'Theta 2', 'X_A', 'Y_A');
for k = idx
    fprintf('%10d %12.2f %12.2f\n', in(k), Ax(k), Ay(k));
end

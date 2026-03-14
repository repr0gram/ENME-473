% computes posture, velocity, and acceleration results
% and generates the plots for Question 4
clc; 
clear;

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
theta1 = deg2rad(180 - atand(70.6/237.2));
alpha = deg2rad(3.2);
beta = deg2rad(3.5);

omega2 = 1;     % rad/s
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
        f = [R2*cos(theta2) + R23*cos(theta23) - R14*cos(theta14) - R1*cos(theta1);
             R2*sin(theta2) + R23*sin(theta23) - R14*sin(theta14) - R1*sin(theta1);
             R2*cos(theta2) + R4*cos(theta23) + R46*cos(theta46) - R36*cos(theta36) - R3*cos(theta14+alpha) - R1*cos(theta1);
             R2*sin(theta2) + R4*sin(theta23) + R46*sin(theta46) - R36*sin(theta36) - R3*sin(theta14+alpha) - R1*sin(theta1);
             R2*cos(theta2) + R4*cos(theta23) + R6*cos(theta46) - R8*cos(theta8) + R7*cos(theta7) - R5*cos(theta36+beta) - R3*cos(theta14+alpha) - R1*cos(theta1);
             R2*sin(theta2) + R4*sin(theta23) + R6*sin(theta46) - R8*sin(theta8) + R7*sin(theta7) - R5*sin(theta36+beta) - R3*sin(theta14+alpha) - R1*sin(theta1)];

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
    b_vel = [R2*sin(theta2);
            -R2*cos(theta2);
             R2*sin(theta2);
            -R2*cos(theta2);
             R2*sin(theta2);
            -R2*cos(theta2)] * omega2;

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
               R2*cos(theta2)*omega2^2 + R4*cos(theta23)*w23^2 + R46*cos(theta46)*w46^2 - R36*cos(theta36)*w36^2 - R3*cos(theta14+alpha)*w14^2;
               R2*sin(theta2)*omega2^2 + R4*sin(theta23)*w23^2 + R46*sin(theta46)*w46^2 - R36*sin(theta36)*w36^2 - R3*sin(theta14+alpha)*w14^2;
               R2*cos(theta2)*omega2^2 + R4*cos(theta23)*w23^2 + R6*cos(theta46)*w46^2 - R8*cos(theta8)*w8^2 + R7*cos(theta7)*w7^2 - R5*cos(theta36+beta)*w36^2 - R3*cos(theta14+alpha)*w14^2;
               R2*sin(theta2)*omega2^2 + R4*sin(theta23)*w23^2 + R6*sin(theta46)*w46^2 - R8*sin(theta8)*w8^2 + R7*sin(theta7)*w7^2 - R5*sin(theta36+beta)*w36^2 - R3*sin(theta14+alpha)*w14^2];

    alpha_vec = J \ rhs_acc;

    alpha3_vals(k) = alpha_vec(2);
    alpha4_vals(k) = alpha_vec(1);
    alpha5_vals(k) = alpha_vec(4);
    alpha6_vals(k) = alpha_vec(3);
    alpha7_vals(k) = alpha_vec(6);
    alpha8_vals(k) = alpha_vec(5);
end

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
title('Angular Velocities vs Input Angle (\omega_2 = 1 rad/s)')
legend('\omega_3','\omega_4','\omega_5','\omega_6','\omega_7','\omega_8','Location','best')
grid on
box on
%exportgraphics(gcf, 'ENME_473_Project_4a_velocity.png', 'Resolution', 600);

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
title('Angular Accelerations vs Input Angle (\omega_2 = 1 rad/s, \alpha_2 = 0)')
legend('\alpha_3','\alpha_4','\alpha_5','\alpha_6','\alpha_7','\alpha_8','Location','best')
grid on
box on
%exportgraphics(gcf, 'ENME_473_Project_4a_acceleration.png', 'Resolution', 600);

%% Write tables to Excel
filename = 'Q4_Results.xlsx';

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

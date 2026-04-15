% You first need to run the Assembly1.slx file found within the same folder
% as this MATLAB file because the output is produced from the .slx file

% ALTERNATIVELY: Drag the out.mat file into your workspace (The out.mat
% file contains results from the Simscape)

% Interim Q4 MATLAB Code Below Labelled as ANALYTICAL

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

t = linspace(0,4.2,121);        % match sim time
theta2_vals = omega2 * t;       % θ = ωt
in = rad2deg(theta2_vals);      % for plotting if needed

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

%% ACTUAL Q5 Comparing ANALYTICAL and SIMULATION values

% --- Extract Simscape data ---
t_sim = out.px.time;

x_sim = out.px.signals.values;
y_sim = out.py.signals.values;

vx_sim = out.vx.signals.values;
vy_sim = out.vy.signals.values;

% Convert to mm
x_sim = x_sim * 1000;
y_sim = y_sim * 1000;
vx_sim = vx_sim * 1000;
vy_sim = vy_sim * 1000;

% Shift Simscape data
x_sim = x_sim - 30;
y_sim = y_sim - 110;

% Limit to 0–4.2 seconds
idx_sim = (t_sim >= 0) & (t_sim <= 4.2);

% (Optional) still compute angle if needed for plotting
theta_sim_deg = rad2deg(t_sim);

theta_sim_deg = theta_sim_deg(idx_sim);
x_sim = x_sim(idx_sim);
y_sim = y_sim(idx_sim);
vx_sim = vx_sim(idx_sim);
vy_sim = vy_sim(idx_sim);

%% POSITION (TRAJECTORY COMPARISON)
figure
plot(Ax, Ay, 'k--', 'LineWidth', 2)   % analytical (dashed black)
hold on
plot(x_sim, y_sim, 'r', 'LineWidth', 2)  % simscape (red)
hold off
xlabel('X Position (mm)')
ylabel('Y Position (mm)')
title('Pin A Trajectory: Analytical vs Simscape')
legend('Analytical','Simscape')
grid on
axis equal

%% VELOCITY COMPARISON (VS ANGLE)
figure
plot(t, Ax_dot, '--', 'Color', [0.6 0 0.6], 'LineWidth', 2)
hold on
plot(t, Ay_dot, '--', 'Color', [0 0.6 0], 'LineWidth', 2)

plot(t_sim(idx_sim), vx_sim, 'Color', [1 0.4 1], 'LineWidth', 2)
plot(t_sim(idx_sim), vy_sim, 'Color', [0.5 1 0], 'LineWidth', 2)

hold off
xlabel('Time (s)')
ylabel('Velocity (mm/s)')
title('Pin A Velocity: Analytical vs Simscape')
legend('Vx Analytical','Vy Analytical','Vx Simscape','Vy Simscape','Location','best')
grid on
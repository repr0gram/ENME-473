% ENME 473 Project - Full Kinematics + Kinetics
% Part I: Position, Velocity, Acceleration, Pin A Path
% Part II: Driving Torque using Instantaneous Power Balance

clc;
clear;
close all;

%  USER INPUTS
omega2 = 1;          % rad/s, constant input angular velocity
alpha2 = 0;          % rad/s^2, constant-speed input assumption

m2 = 0.3474;         % kg, input link mass  (CHANGE to your team's value)
m8 = 1.154112;       % kg, output link mass (CHANGE to your team's value)
g  = 9.81;           % m/s^2

%%  KNOWN GEOMETRY  (all lengths in mm)
R1  = sqrt(237.2^2 + 70.6^2);
R2  = 272.4;
R23 = 236.9;
R14 = 279.0;
R3  = 489.9;
R4  = 539.8;
R36 = 342.1;
R46 = 179.3;
R5  = 595.5;
R6  = 378.9;
R7  = 198.2;
R8  = 254.4;
RA  = 336.9;               % distance from joint 7/8 to Pin A along link 8

theta1 = deg2rad(180 - atand(70.6/237.2));
alpha  = deg2rad(3.2);
beta   = deg2rad(3.5);
gamma  = deg2rad(1.8);

RAtot = R8 + RA;           % total distance from link-8 pivot to Pin A (mm)

% Mass moments of inertia about CG  (kg·m²)
R2_m = R2/1000;            % link 2 length in metres
R8_m = R8/1000;            % link 8 length in metres
IG2  = m2 * R2_m^2 / 12;
IG8  = m8 * R8_m^2 / 12;

%% =========================
%  INPUT RANGE
in            = 0:120;              % degrees, 1-degree increments
theta2_vals   = deg2rad(in);        % radians
N             = length(theta2_vals);

%% =========================
%  STORAGE ARRAYS
% --- solved angles (rad) ---
theta23_raw = zeros(1,N);
theta14_raw = zeros(1,N);
theta46_raw = zeros(1,N);
theta36_raw = zeros(1,N);
theta8_raw  = zeros(1,N);
theta7_raw  = zeros(1,N);

% --- angular velocities (rad/s) ---
omega23_vals = zeros(1,N);
omega14_vals = zeros(1,N);
omega46_vals = zeros(1,N);
omega36_vals = zeros(1,N);
omega8_vals  = zeros(1,N);
omega7_vals  = zeros(1,N);

% --- angular accelerations (rad/s²) ---
alpha23_vals = zeros(1,N);
alpha14_vals = zeros(1,N);
alpha46_vals = zeros(1,N);
alpha36_vals = zeros(1,N);
alpha8_vals  = zeros(1,N);
alpha7_vals  = zeros(1,N);

% --- Pin A position (mm), velocity (mm/s), acceleration (mm/s²) ---
Ax      = zeros(1,N);
Ay      = zeros(1,N);
Ax_dot  = zeros(1,N);
Ay_dot  = zeros(1,N);
Ax_ddot = zeros(1,N);
Ay_ddot = zeros(1,N);

% --- link 2 CG kinematics (mm, mm/s, mm/s²) ---
xG2  = zeros(1,N);
yG2  = zeros(1,N);
vxG2 = zeros(1,N);
vyG2 = zeros(1,N);
axG2 = zeros(1,N);
ayG2 = zeros(1,N);

% --- link 8 CG kinematics (mm, mm/s, mm/s²) ---
xG8  = zeros(1,N);
yG8  = zeros(1,N);
vxG8 = zeros(1,N);
vyG8 = zeros(1,N);
axG8 = zeros(1,N);
ayG8 = zeros(1,N);

% --- driving torque ---
T2_drive = zeros(1,N);

%% =========================
%  NEWTON-RAPHSON SETTINGS
%  =========================
% Initial guess chosen to converge to the correct assembly / branch
% for the desired Pin A path.
x = deg2rad([167.5; -4.0; 10.7; 165.4; -4.8; -171.3]);

tol_f   = 1e-8;
tol_x   = 1e-8;
maxIter = 100;

%% =========================
%  MAIN LOOP
for k = 1:N

    theta2 = theta2_vals(k);

    % -------------------------------------------------------
    % (1) POSITION ANALYSIS  —  Newton-Raphson
    for iter = 1:maxIter

        theta23 = x(1);
        theta14 = x(2);
        theta46 = x(3);
        theta36 = x(4);
        theta8  = x(5);
        theta7  = x(6);

        % Vector-loop residuals
        f = [ R2*cos(theta2) + R23*cos(theta23) - R14*cos(theta14) - R1*cos(theta1);
              R2*sin(theta2) + R23*sin(theta23) - R14*sin(theta14) - R1*sin(theta1);
              R2*cos(theta2) + R4*cos(theta23+gamma)  + R46*cos(theta46) - R36*cos(theta36) - R3*cos(theta14+alpha) - R1*cos(theta1);
              R2*sin(theta2) + R4*sin(theta23+gamma)  + R46*sin(theta46) - R36*sin(theta36) - R3*sin(theta14+alpha) - R1*sin(theta1);
              R2*cos(theta2) + R4*cos(theta23+gamma)  + R6*cos(theta46)  - R8*cos(theta8)   + R7*cos(theta7) - R5*cos(theta36+beta) - R3*cos(theta14+alpha) - R1*cos(theta1);
              R2*sin(theta2) + R4*sin(theta23+gamma)  + R6*sin(theta46)  - R8*sin(theta8)   + R7*sin(theta7) - R5*sin(theta36+beta) - R3*sin(theta14+alpha) - R1*sin(theta1)];

        % Jacobian  ∂f/∂[theta23, theta14, theta46, theta36, theta8, theta7]
        J = [-R23*sin(theta23),  R14*sin(theta14),             0,                   0,              0,             0;
              R23*cos(theta23), -R14*cos(theta14),             0,                   0,              0,             0;
             -R4*sin(theta23+gamma),   R3*sin(theta14+alpha),  -R46*sin(theta46),   R36*sin(theta36),     0,             0;
              R4*cos(theta23+gamma),  -R3*cos(theta14+alpha),   R46*cos(theta46),  -R36*cos(theta36),     0,             0;
             -R4*sin(theta23+gamma),   R3*sin(theta14+alpha),   -R6*sin(theta46),   R5*sin(theta36+beta),  R8*sin(theta8), -R7*sin(theta7);
              R4*cos(theta23+gamma),  -R3*cos(theta14+alpha),    R6*cos(theta46),  -R5*cos(theta36+beta), -R8*cos(theta8),  R7*cos(theta7)];

        dx    = J \ f;
        x_new = x - dx;

        % Wrap angles to [-π, π] to stay in the correct branch
        x_new = atan2(sin(x_new), cos(x_new));

        if norm(f, inf) < tol_f && norm(x_new - x, inf) < tol_x
            x = x_new;
            break;
        end

        x = x_new;
    end

    if iter == maxIter
        error('Newton-Raphson did not converge at theta2 = %g deg.', in(k));
    end

    theta23 = x(1);
    theta14 = x(2);
    theta46 = x(3);
    theta36 = x(4);
    theta8  = x(5);
    theta7  = x(6);

    % Store raw posture results
    theta23_raw(k) = theta23;
    theta14_raw(k) = theta14;
    theta46_raw(k) = theta46;
    theta36_raw(k) = theta36;
    theta8_raw(k)  = theta8;
    theta7_raw(k)  = theta7;

    % Rebuild Jacobian at converged posture (reused for velocity & acceleration)
    J = [-R23*sin(theta23),  R14*sin(theta14),             0,                   0,              0,             0;
          R23*cos(theta23), -R14*cos(theta14),             0,                   0,              0,             0;
         -R4*sin(theta23+gamma),   R3*sin(theta14+alpha),  -R46*sin(theta46),   R36*sin(theta36),     0,             0;
          R4*cos(theta23+gamma),  -R3*cos(theta14+alpha),   R46*cos(theta46),  -R36*cos(theta36),     0,             0;
         -R4*sin(theta23+gamma),   R3*sin(theta14+alpha),   -R6*sin(theta46),   R5*sin(theta36+beta),  R8*sin(theta8), -R7*sin(theta7);
          R4*cos(theta23+gamma),  -R3*cos(theta14+alpha),    R6*cos(theta46),  -R5*cos(theta36+beta), -R8*cos(theta8),  R7*cos(theta7)];

    % -------------------------------------------------------
    % (2) VELOCITY ANALYSIS   J * ω_vec = b_vel
    %     b_vel comes from d/dt of each position loop equation,
    %     isolating the omega2 contribution to the RHS.
    % -------------------------------------------------------
    b_vel = omega2 * [ R2*sin(theta2);
                      -R2*cos(theta2);
                       R2*sin(theta2);
                      -R2*cos(theta2);
                       R2*sin(theta2);
                      -R2*cos(theta2)];

    omega_vec = J \ b_vel;
    % omega_vec = [omega23; omega14; omega46; omega36; omega8; omega7]

    omega23 = omega_vec(1);
    omega14 = omega_vec(2);
    omega46 = omega_vec(3);
    omega36 = omega_vec(4);
    omega8  = omega_vec(5);
    omega7  = omega_vec(6);

    omega23_vals(k) = omega23;
    omega14_vals(k) = omega14;
    omega46_vals(k) = omega46;
    omega36_vals(k) = omega36;
    omega8_vals(k)  = omega8;
    omega7_vals(k)  = omega7;

    % -------------------------------------------------------
    % (3) ACCELERATION ANALYSIS   J * α_vec = rhs_acc
    %     RHS = centripetal terms (all −R·ω²·[cos; sin] contributions
    %     moved to the right-hand side after differentiating velocities).
    rhs_acc = [ R2*cos(theta2)*omega2^2  + R23*cos(theta23)*omega23^2  - R14*cos(theta14)*omega14^2;
                R2*sin(theta2)*omega2^2  + R23*sin(theta23)*omega23^2  - R14*sin(theta14)*omega14^2;
                R2*cos(theta2)*omega2^2  + R4*cos(theta23+gamma)*omega23^2   + R46*cos(theta46)*omega46^2 - R36*cos(theta36)*omega36^2 - R3*cos(theta14+alpha)*omega14^2;
                R2*sin(theta2)*omega2^2  + R4*sin(theta23+gamma)*omega23^2   + R46*sin(theta46)*omega46^2 - R36*sin(theta36)*omega36^2 - R3*sin(theta14+alpha)*omega14^2;
                R2*cos(theta2)*omega2^2  + R4*cos(theta23+gamma)*omega23^2   + R6*cos(theta46)*omega46^2  - R8*cos(theta8)*omega8^2   + R7*cos(theta7)*omega7^2 - R5*cos(theta36+beta)*omega36^2 - R3*cos(theta14+alpha)*omega14^2;
                R2*sin(theta2)*omega2^2  + R4*sin(theta23+gamma)*omega23^2   + R6*sin(theta46)*omega46^2  - R8*sin(theta8)*omega8^2   + R7*sin(theta7)*omega7^2 - R5*sin(theta36+beta)*omega36^2 - R3*sin(theta14+alpha)*omega14^2];

    alpha_vec = J \ rhs_acc;
    % alpha_vec = [alpha23; alpha14; alpha46; alpha36; alpha8; alpha7]

    alpha23 = alpha_vec(1);
    alpha14 = alpha_vec(2);
    alpha46 = alpha_vec(3);
    alpha36 = alpha_vec(4);
    alpha8  = alpha_vec(5);
    alpha7  = alpha_vec(6);

    alpha23_vals(k) = alpha23;
    alpha14_vals(k) = alpha14;
    alpha46_vals(k) = alpha46;
    alpha36_vals(k) = alpha36;
    alpha8_vals(k)  = alpha8;
    alpha7_vals(k)  = alpha7;

    % -------------------------------------------------------
    % (4) PIN A  POSITION, VELOCITY, ACCELERATION
    Ax(k) = R2*cos(theta2) + R4*cos(theta23+gamma) + R6*cos(theta46) - RAtot*cos(theta8);
    Ay(k) = R2*sin(theta2) + R4*sin(theta23+gamma) + R6*sin(theta46) - RAtot*sin(theta8);

    Ax_dot(k) = -R2*sin(theta2)*omega2 ...
                -R4*sin(theta23+gamma)*omega23 ...
                -R6*sin(theta46)*omega46 ...
                +RAtot*sin(theta8)*omega8;

    Ay_dot(k) =  R2*cos(theta2)*omega2 ...
                +R4*cos(theta23+gamma)*omega23 ...
                +R6*cos(theta46)*omega46 ...
                -RAtot*cos(theta8)*omega8;

    Ax_ddot(k) = -R2*cos(theta2)*omega2^2 ...
                 -R4*cos(theta23+gamma)*omega23^2 - R4*sin(theta23+gamma)*alpha23 ...
                 -R6*cos(theta46)*omega46^2 - R6*sin(theta46)*alpha46 ...
                 +RAtot*cos(theta8)*omega8^2 + RAtot*sin(theta8)*alpha8;

    % CORRECTED: was "+RAtot*sin(theta8)*alpha8" — now "-RAtot*cos(theta8)*alpha8"
    Ay_ddot(k) = -R2*sin(theta2)*omega2^2 ...
                 -R4*sin(theta23+gamma)*omega23^2 + R4*cos(theta23+gamma)*alpha23 ...
                 -R6*sin(theta46)*omega46^2 + R6*cos(theta46)*alpha46 ...
                 +RAtot*sin(theta8)*omega8^2 - RAtot*cos(theta8)*alpha8;

    % -------------------------------------------------------
    % (5) LINK 2 CG KINEMATICS
    %     CG at midpoint of link 2 (from ground pivot O2):
    %       xG2 = (R2/2)·cos(θ2),  yG2 = (R2/2)·sin(θ2)
    xG2(k)  =  (R2/2)*cos(theta2);
    yG2(k)  =  (R2/2)*sin(theta2);

    vxG2(k) = -(R2/2)*sin(theta2)*omega2;
    vyG2(k)  =  (R2/2)*cos(theta2)*omega2;

    axG2(k) = -(R2/2)*cos(theta2)*omega2^2 - (R2/2)*sin(theta2)*alpha2;
    ayG2(k) = -(R2/2)*sin(theta2)*omega2^2 + (R2/2)*cos(theta2)*alpha2;

    % -------------------------------------------------------
    % (6) LINK 8 CG KINEMATICS

    xP = R2*cos(theta2) + R4*cos(theta23+gamma) + R6*cos(theta46);
    yP = R2*sin(theta2) + R4*sin(theta23+gamma) + R6*sin(theta46);

    vxP = -R2*sin(theta2)*omega2 ...
          -R4*sin(theta23+gamma)*omega23 ...
          -R6*sin(theta46)*omega46;

    vyP =  R2*cos(theta2)*omega2 ...
          +R4*cos(theta23+gamma)*omega23 ...
          +R6*cos(theta46)*omega46;

    axP = -R2*cos(theta2)*omega2^2 ...
          -R2*sin(theta2)*alpha2 ...
          -R4*cos(theta23+gamma)*omega23^2 - R4*sin(theta23+gamma)*alpha23 ...
          -R6*cos(theta46)*omega46^2 - R6*sin(theta46)*alpha46;

    ayP = -R2*sin(theta2)*omega2^2 ...
          +R2*cos(theta2)*alpha2 ...
          -R4*sin(theta23+gamma)*omega23^2 + R4*cos(theta23+gamma)*alpha23 ...
          -R6*sin(theta46)*omega46^2 + R6*cos(theta46)*alpha46;

    xG8(k)  = xP - (R8/2)*cos(theta8);
    yG8(k)  = yP - (R8/2)*sin(theta8);

    vxG8(k) = vxP + (R8/2)*sin(theta8)*omega8;
    vyG8(k) = vyP - (R8/2)*cos(theta8)*omega8;

    axG8(k) = axP + (R8/2)*cos(theta8)*omega8^2 + (R8/2)*sin(theta8)*alpha8;
    ayG8(k) = ayP + (R8/2)*sin(theta8)*omega8^2 - (R8/2)*cos(theta8)*alpha8;

    % -------------------------------------------------------
    % (7) DRIVING TORQUE  —  Instantaneous Power Balance
    % Convert mm-based quantities to SI (m, m/s, m/s²)
    vG2x_si = vxG2(k)/1000;   vG2y_si = vyG2(k)/1000;
    aG2x_si = axG2(k)/1000;   aG2y_si = ayG2(k)/1000;

    vG8x_si = vxG8(k)/1000;   vG8y_si = vyG8(k)/1000;
    aG8x_si = axG8(k)/1000;   aG8y_si = ayG8(k)/1000;

    % Rate of change of kinetic energy for each link  (W)
    Kdot2 = m2*(aG2x_si*vG2x_si + aG2y_si*vG2y_si) + IG2*alpha2*omega2;
    Kdot8 = m8*(aG8x_si*vG8x_si + aG8y_si*vG8y_si) + IG8*alpha8*omega8;

    % T2·ω2 = Kdot_total + m·g·ẏG  (gravity power — positive when CG rises)
    T2_drive(k) = (Kdot2 + Kdot8 ...
                   + m2*g*vG2y_si ...
                   + m8*g*vG8y_si) / omega2;
end

%% =========================
%  POST-PROCESS ANGLES
theta23_deg = rad2deg(unwrap(theta23_raw));
theta14_deg = rad2deg(unwrap(theta14_raw));
theta46_deg = rad2deg(unwrap(theta46_raw));
theta36_deg = rad2deg(unwrap(theta36_raw));
theta8_deg  = rad2deg(unwrap(theta8_raw));
theta7_deg  = rad2deg(unwrap(theta7_raw));

% Physical moving-link angles for reporting
theta3_deg      = theta14_deg + rad2deg(alpha);
theta4_deg      = theta23_deg + rad2deg(gamma);
theta5_deg      = theta36_deg + rad2deg(beta);
theta6_deg      = theta46_deg;
theta7_plot_deg = theta7_deg  + 180;
theta8_plot_deg = theta8_deg  + 180;

% Physical link angular velocities
omega3_vals      = omega14_vals;
omega4_vals      = omega23_vals;
omega5_vals      = omega36_vals;
omega6_vals      = omega46_vals;
omega7_plot_vals = omega7_vals;
omega8_plot_vals = omega8_vals;

% Physical link angular accelerations
alpha3_plot_vals = alpha14_vals;
alpha4_plot_vals = alpha23_vals;
alpha5_plot_vals = alpha36_vals;
alpha6_plot_vals = alpha46_vals;
alpha7_plot_vals = alpha7_vals;
alpha8_plot_vals = alpha8_vals;

%% =========================
%  PLOTS
scriptDir = fileparts(mfilename('fullpath'));

% 1) Link angles vs input angle
figure
plot(in, theta3_deg, 'LineWidth', 2); hold on
plot(in, theta4_deg, 'LineWidth', 2)
plot(in, theta5_deg, 'LineWidth', 2)
plot(in, theta6_deg, 'LineWidth', 2)
plot(in, theta7_plot_deg, 'LineWidth', 2)
plot(in, theta8_plot_deg, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Link Angle (deg)')
title('Link Angles vs Input Angle')
legend('\theta_3','\theta_4','\theta_5','\theta_6','\theta_7','\theta_8','Location','best')
grid on; box on
exportgraphics(gcf, fullfile(scriptDir, 'ENME473_LinkAngles.png'), 'Resolution', 600);

% 2) Path traced by Pin A
figure
plot(Ax, Ay, 'b-', 'LineWidth', 2); hold on
plot(Ax(1),   Ay(1),   'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', ...
    'DisplayName', '\theta_2 = 0°')
plot(Ax(end), Ay(end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
    'DisplayName', '\theta_2 = 120°')
markAngles = [30 60 90];
for n = 1:length(markAngles)
    idx = markAngles(n) + 1;
    plot(Ax(idx), Ay(idx), 'kd', 'MarkerSize', 8, 'MarkerFaceColor', 'k', ...
        'HandleVisibility', 'off')
    text(Ax(idx)+10, Ay(idx)+10, sprintf('\\theta_2 = %d°', markAngles(n)), ...
        'FontSize', 11)
end
xlabel('X (mm)'); ylabel('Y (mm)')
title('Path Traced by Pin A')
legend('Path of Pin A', '\theta_2 = 0°', '\theta_2 = 120°', 'Location', 'best')
grid on; box on; axis equal
exportgraphics(gcf, fullfile(scriptDir, 'ENME473_PinA_Path.png'), 'Resolution', 600);

% 3) Angular velocities vs input angle
figure
plot(in, omega3_vals, 'LineWidth', 2); hold on
plot(in, omega4_vals, 'LineWidth', 2)
plot(in, omega5_vals, 'LineWidth', 2)
plot(in, omega6_vals, 'LineWidth', 2)
plot(in, omega7_plot_vals, 'LineWidth', 2)
plot(in, omega8_plot_vals, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Angular Velocity (rad/s)')
title('Angular Velocities vs Input Angle')
legend('\omega_3','\omega_4','\omega_5','\omega_6','\omega_7','\omega_8','Location','best')
grid on; box on
exportgraphics(gcf, fullfile(scriptDir, 'ENME473_AngularVelocities.png'), 'Resolution', 600);

% 4) Angular accelerations vs input angle
figure
plot(in, alpha3_plot_vals, 'LineWidth', 2); hold on
plot(in, alpha4_plot_vals, 'LineWidth', 2)
plot(in, alpha5_plot_vals, 'LineWidth', 2)
plot(in, alpha6_plot_vals, 'LineWidth', 2)
plot(in, alpha7_plot_vals, 'LineWidth', 2)
plot(in, alpha8_plot_vals, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Angular Acceleration (rad/s^2)')
title('Angular Accelerations vs Input Angle')
legend('\alpha_3','\alpha_4','\alpha_5','\alpha_6','\alpha_7','\alpha_8','Location','best')
grid on; box on
exportgraphics(gcf, fullfile(scriptDir, 'ENME473_AngularAccelerations.png'), 'Resolution', 600);

% 5) Pin A velocity vs input angle
figure
plot(in, Ax_dot, 'LineWidth', 2); hold on
plot(in, Ay_dot, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Velocity (mm/s)')
title('Pin A Velocity vs Input Angle')
legend('v_{Ax}', 'v_{Ay}', 'Location', 'best')
grid on; box on
exportgraphics(gcf, fullfile(scriptDir, 'ENME473_PinA_Velocity.png'), 'Resolution', 600);

% 6) Pin A acceleration vs input angle
figure
plot(in, Ax_ddot, 'LineWidth', 2); hold on
plot(in, Ay_ddot, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Acceleration (mm/s^2)')
title('Pin A Acceleration vs Input Angle')
legend('a_{Ax}', 'a_{Ay}', 'Location', 'best')
grid on; box on
exportgraphics(gcf, fullfile(scriptDir, 'ENME473_PinA_Acceleration.png'), 'Resolution', 600);

% 7) Driving torque vs input angle
figure
plot(in, T2_drive, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Driving Torque T_2 (N\cdotm)')
title(sprintf('Driving Torque vs Input Angle  (m_2 = %.4f kg,  m_8 = %.4f kg)', m2, m8))
grid on; box on
exportgraphics(gcf, fullfile(scriptDir, 'ENME473_DrivingTorque.png'), 'Resolution', 600);

%% =========================
%  WRITE RESULTS TO EXCEL
filename = fullfile(scriptDir, 'ENME473_FullResults.xlsx');

T_angles = table(in', theta3_deg', theta4_deg', theta5_deg', theta6_deg', ...
    theta7_plot_deg', theta8_plot_deg', ...
    'VariableNames', {'Theta2_deg','Theta3_deg','Theta4_deg','Theta5_deg', ...
                      'Theta6_deg','Theta7_deg','Theta8_deg'});
writetable(T_angles, filename, 'Sheet', 'Link Angles');

T_vel = table(in', omega3_vals', omega4_vals', omega5_vals', omega6_vals', ...
    omega7_plot_vals', omega8_plot_vals', ...
    'VariableNames', {'Theta2_deg','omega3','omega4','omega5','omega6','omega7','omega8'});
writetable(T_vel, filename, 'Sheet', 'Angular Velocities');

T_acc = table(in', alpha3_plot_vals', alpha4_plot_vals', alpha5_plot_vals', ...
    alpha6_plot_vals', alpha7_plot_vals', alpha8_plot_vals', ...
    'VariableNames', {'Theta2_deg','alpha3','alpha4','alpha5','alpha6','alpha7','alpha8'});
writetable(T_acc, filename, 'Sheet', 'Angular Accelerations');

T_pinA = table(in', Ax', Ay', Ax_dot', Ay_dot', Ax_ddot', Ay_ddot', ...
    'VariableNames', {'Theta2_deg','X_A_mm','Y_A_mm', ...
                      'Vx_A_mm_s','Vy_A_mm_s','Ax_A_mm_s2','Ay_A_mm_s2'});
writetable(T_pinA, filename, 'Sheet', 'Pin A');

T_cg = table(in', xG2', yG2', vxG2', vyG2', axG2', ayG2', ...
    xG8', yG8', vxG8', vyG8', axG8', ayG8', ...
    'VariableNames', {'Theta2_deg', ...
    'xG2_mm','yG2_mm','vxG2_mm_s','vyG2_mm_s','axG2_mm_s2','ayG2_mm_s2', ...
    'xG8_mm','yG8_mm','vxG8_mm_s','vyG8_mm_s','axG8_mm_s2','ayG8_mm_s2'});
writetable(T_cg, filename, 'Sheet', 'CG Kinematics');

T_torque = table(in', T2_drive', ...
    'VariableNames', {'Theta2_deg','DrivingTorque_Nm'});
writetable(T_torque, filename, 'Sheet', 'Driving Torque');

%% =========================
%  DEDICATED TORQUE EXCEL FILE
%  Exports theta2 vs T2 to a standalone, clearly labelled spreadsheet.
torque_filename = fullfile(scriptDir, 'ENME473_DrivingTorque.xlsx');

% --- Main data table: theta2 and T2 for every 1-degree increment ----------
T_torque_export = table(in', T2_drive', ...
    'VariableNames', {'Theta2_deg', 'DrivingTorque_Nm'});
writetable(T_torque_export, torque_filename, 'Sheet', 'Torque vs Theta');

% --- Summary statistics appended below the data via writecell -------------
[maxT, idxMaxT] = max(T2_drive);
[minT, idxMinT] = min(T2_drive);

summaryData = {
    '',           '';
    'Summary',    '';
    'omega2 (rad/s)',            omega2;
    'alpha2 (rad/s^2)',          alpha2;
    'Input link mass m2 (kg)',   m2;
    'Output link mass m8 (kg)',  m8;
    'Max T2 (N·m)',              maxT;
    'Max T2 at theta2 (deg)',    in(idxMaxT);
    'Min T2 (N·m)',              minT;
    'Min T2 at theta2 (deg)',    in(idxMinT);
};
% Write summary below the data (data occupies rows 1 to N+1, so start at N+3)
startRow = length(in) + 3;
writecell(summaryData, torque_filename, 'Sheet', 'Torque vs Theta', ...
    'Range', sprintf('A%d', startRow));

fprintf('\nDriving torque data written to: %s\n', torque_filename);
fprintf('  Sheet : Torque vs Theta\n');
fprintf('  Rows  : %d data rows (theta2 = 0 to %d deg, 1-deg steps)\n', ...
    length(in), in(end));
fprintf('  Max T2 = %.6f N·m at theta2 = %d deg\n', maxT, in(idxMaxT));
fprintf('  Min T2 = %.6f N·m at theta2 = %d deg\n', minT, in(idxMinT));

%% =========================
%  PRINT KEY RESULTS
idx = [1 31 61 91 121];   % theta2 = 0, 30, 60, 90, 120 deg

fprintf('\n====================================================\n');
fprintf('ENME 473 PROJECT RESULTS SUMMARY\n');
fprintf('====================================================\n');
fprintf('Input angular velocity  omega2 = %.4f rad/s\n', omega2);
fprintf('Input angular accel.   alpha2 = %.4f rad/s^2\n', alpha2);
fprintf('Input  link mass           m2 = %.4f kg\n', m2);
fprintf('Output link mass           m8 = %.4f kg\n', m8);
fprintf('Results written to %s\n', filename);

fprintf('\n=== Pin A Position (mm) ===\n');
fprintf('%10s %14s %14s\n', 'Theta2', 'X_A', 'Y_A');
for k = idx
    fprintf('%10d %14.4f %14.4f\n', in(k), Ax(k), Ay(k));
end

fprintf('\n=== Pin A Velocity (mm/s) ===\n');
fprintf('%10s %14s %14s\n', 'Theta2', 'Vx_A', 'Vy_A');
for k = idx
    fprintf('%10d %14.4f %14.4f\n', in(k), Ax_dot(k), Ay_dot(k));
end

fprintf('\n=== Pin A Acceleration (mm/s^2) ===\n');
fprintf('%10s %14s %14s\n', 'Theta2', 'Ax_A', 'Ay_A');
for k = idx
    fprintf('%10d %14.4f %14.4f\n', in(k), Ax_ddot(k), Ay_ddot(k));
end

fprintf('\n=== Driving Torque (N·m) ===\n');
fprintf('%10s %14s\n', 'Theta2', 'T2');
for k = idx
    fprintf('%10d %14.6f\n', in(k), T2_drive(k));
end

[maxTorque, idxMax] = max(T2_drive);
[minTorque, idxMin] = min(T2_drive);

fprintf('\nMaximum driving torque = %.6f N·m at theta2 = %d deg\n', maxTorque, in(idxMax));
fprintf('Minimum driving torque = %.6f N·m at theta2 = %d deg\n', minTorque, in(idxMin));
fprintf('Negative torque means the actuator is braking/resisting the motion.\n');
% You first need to run the Assembly1.slx file found within the same folder
% as this MATLAB file because the output is produced from the .slx file

% Extract data from Simulation Output
t = out.px.time;

x = out.px.signals.values;
y = out.py.signals.values;

vx_data = out.vx.signals.values;
vy_data = out.vy.signals.values;

% Limit to t = [0, 10]
idx = (t >= 0) & (t <= 10);

t = t(idx);
x = x(idx);
y = y(idx);
vx_data = vx_data(idx);
vy_data = vy_data(idx);

% --- Trajectory (Y vs X) ---
figure
plot(x, y, 'r', 'LineWidth', 2)   % red line, thicker
xlabel('X Position (m)')
ylabel('Y Position (m)')
title('Position (Y vs X)')
grid on
axis equal

% --- Velocity ---
figure
plot(t, vx_data, 'Color', [1 0.4 1], 'LineWidth', 2)  % pink line
hold on
plot(t, vy_data, 'Color', [0.5 1 0], 'LineWidth', 2)  % lime green line
hold off
xlabel('Time (s)')
ylabel('Velocity (m/s)')
legend('Vx','Vy')
title('Velocity vs Time')
grid on
% ENME 473 Project - Part III: Design (Question 7)
% finds the "goal" Pin A position at theta2 = 120 deg,
% computes the enclosing box at theta2 = 0 deg,
% and systematically varies design parameters to minimize box area.

clc; clear; close all;
scriptDir = fileparts(mfilename('fullpath'));

%% original geometry (mm)
R1_orig  = sqrt(237.2^2 + 70.6^2);
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
R8_orig  = 254.4;
RA_orig  = 336.9;

theta1_orig = deg2rad(180 - atand(70.6/237.2));
alpha  = deg2rad(3.2);
beta   = deg2rad(3.5);
gamma  = deg2rad(1.8);

% ground pivot O4 original position (mm)
xO4_orig = -237.2;
yO4_orig = 70.6;

%% solve original design
fprintf("=== ORIGINAL DESIGN ===\n");

% initial guess for NR (converges to correct branch)
x0 = deg2rad([167.5; -4.0; 10.7; 165.4; -4.8; -171.3]);

% solve full range 0-120 deg
[allAngles, pinA, boxInfo] = solveDesign(R2, R23, R14, R3, R4, R36, R46, ...
    R5, R6, R7, R8_orig, RA_orig, R1_orig, theta1_orig, alpha, beta, gamma, x0);

% goal: Pin A at theta2 = 120 deg (last entry)
goalX = pinA.Ax(end);
goalY = pinA.Ay(end);
fprintf("Pin A at theta2 = 120 deg (GOAL): (%.2f, %.2f) mm\n", goalX, goalY);

% enclosing box at theta2 = 0 deg
fprintf("Enclosing box at theta2 = 0 deg:\n");
fprintf("  Width  = %.2f mm\n", boxInfo.width);
fprintf("  Height = %.2f mm\n", boxInfo.height);
fprintf("  Area   = %.2f mm^2\n", boxInfo.area);
fprintf("  X range: [%.2f, %.2f]\n", boxInfo.xmin, boxInfo.xmax);
fprintf("  Y range: [%.2f, %.2f]\n", boxInfo.ymin, boxInfo.ymax);

origArea = boxInfo.area;

%% plot original Pin A path and mechanism at theta2 = 0
figure
subplot(1, 2, 1)
plotMechanism(boxInfo.joints, pinA.Ax(1), pinA.Ay(1), boxInfo);
title("Mechanism at \theta_2 = 0° (Original)")

subplot(1, 2, 2)
plot(pinA.Ax, pinA.Ay, "b-", LineWidth=2)
hold on
plot(pinA.Ax(1), pinA.Ay(1), "go", MarkerSize=10, MarkerFaceColor="g")
plot(goalX, goalY, "rs", MarkerSize=10, MarkerFaceColor="r")
xlabel("X (mm)"); ylabel("Y (mm)")
title("Pin A Path (Original)")
legend("Path", "\theta_2 = 0°", "\theta_2 = 120° (Goal)", Location="best")
grid on; axis equal; box on
exportgraphics(gcf, fullfile(scriptDir, "ENME473_Q7_Original.png"), Resolution=600);

%% plot original mechanism and path every 15 deg from 0 to 120
plotAnglesFigure(boxInfo.allJoints, pinA, 0:15:120, goalX, goalY, "Original");
exportgraphics(gcf, fullfile(scriptDir, "ENME473_Q7_Original_Angles.png"), Resolution=600);

%% summary table: box dimensions and area every 5 deg (original)
xlsxPath = fullfile(scriptDir, "ENME473_Q7_BoxDimensions.xlsx");
if isfile(xlsxPath)
    delete(xlsxPath) % start fresh each run
end
printAngleTable(boxInfo.allJoints, "Original", xlsxPath);

%% design parameter search
% parameters to vary:
%   (a) ground pivot O4 location: (xO4, yO4)
%   (b) RA: distance from Pin A to joint 7/8
%   (c) R8: distance from joint 6/8 to joint 7/8
%
% strategy: grid search around original values, then report best designs.
% for each candidate, check if Pin A passes near the "goal" at any theta2.

goalTol = 15; % mm tolerance for reaching goal position

% search ranges (vary +-20% around original, coarse then fine)
xO4_range = linspace(xO4_orig - 50, xO4_orig + 50, 11);
yO4_range = linspace(yO4_orig - 30, yO4_orig + 30, 7);
RA_range  = linspace(RA_orig*0.7, RA_orig*1.1, 9);
R8_range  = linspace(R8_orig*0.8, R8_orig*1.1, 7);

% get the converged angles from the original design to use as adaptive guess
x0_base = allAngles.raw(:, 1); % converged solution at theta2 = 0 for original

% storage for results
results = [];
nTotal = length(xO4_range) * length(yO4_range) * length(RA_range) * length(R8_range);
fprintf("\n=== DESIGN SEARCH ===\n");
fprintf("Evaluating %d design candidates...\n", nTotal);

count = 0;
for iX = 1:length(xO4_range)
    for iY = 1:length(yO4_range)
        xO4 = xO4_range(iX);
        yO4 = yO4_range(iY);

        R1_new = sqrt(xO4^2 + yO4^2);
        theta1_new = atan2(yO4, xO4);

        for iRA = 1:length(RA_range)
            for iR8 = 1:length(R8_range)
                RA_new = RA_range(iRA);
                R8_new = R8_range(iR8);
                count = count + 1;

                try
                    % use converged original solution as initial guess
                    % so NR can track the solution as parameters shift
                    [~, pA, bx] = solveDesign(R2, R23, R14, R3, R4, R36, R46, ...
                        R5, R6, R7, R8_new, RA_new, R1_new, theta1_new, ...
                        alpha, beta, gamma, x0_base);

                    % check if Pin A passes near goal at any theta2
                    dist = sqrt((pA.Ax - goalX).^2 + (pA.Ay - goalY).^2);
                    [minDist, idxMin] = min(dist);

                    if minDist < goalTol && bx.area < origArea
                        r.xO4    = xO4;
                        r.yO4    = yO4;
                        r.RA     = RA_new;
                        r.R8     = R8_new;
                        r.area   = bx.area;
                        r.width  = bx.width;
                        r.height = bx.height;
                        r.goalDist = minDist;
                        r.goalTheta2 = idxMin - 1; % degrees
                        r.goalAx = pA.Ax(idxMin);
                        r.goalAy = pA.Ay(idxMin);
                        results = [results; r]; %#ok<AGROW>
                    end
                catch
                    % NR did not converge for this design — skip
                end
            end
        end
    end
end

fprintf("Search complete. Found %d designs with smaller box reaching the goal.\n", length(results));

%% report results
if ~isempty(results)
    % sort by area
    areas = [results.area];
    [~, sortIdx] = sort(areas);
    results = results(sortIdx);

    % report best design
    best = results(1);
    fprintf("\n=== BEST DESIGN ===\n");
    fprintf("Ground pivot O4: (%.2f, %.2f) mm\n", best.xO4, best.yO4);
    fprintf("RA = %.2f mm,  R8 = %.2f mm\n", best.RA, best.R8);
    fprintf("Enclosing box: %.2f x %.2f mm,  Area = %.2f mm^2\n", ...
        best.width, best.height, best.area);
    fprintf("Area reduction: %.2f%%\n", (1 - best.area/origArea)*100);
    fprintf("Goal reached at theta2 = %d deg (distance = %.2f mm)\n", ...
        best.goalTheta2, best.goalDist);
    fprintf("Pin A at goal: (%.2f, %.2f) mm\n", best.goalAx, best.goalAy);

    % plot best design
    [~, pA_best, bx_best] = solveDesign(R2, R23, R14, R3, R4, R36, R46, ...
        R5, R6, R7, best.R8, best.RA, ...
        sqrt(best.xO4^2 + best.yO4^2), atan2(best.yO4, best.xO4), ...
        alpha, beta, gamma, x0_base);

    figure
    subplot(1, 2, 1)
    plotMechanism(bx_best.joints, pA_best.Ax(1), pA_best.Ay(1), bx_best);
    title("Mechanism at \theta_2 = 0° (Best Design)")

    subplot(1, 2, 2)
    plot(pA_best.Ax, pA_best.Ay, "b-", LineWidth=2)
    hold on
    plot(pA_best.Ax(1), pA_best.Ay(1), "go", MarkerSize=10, MarkerFaceColor="g")
    plot(goalX, goalY, "rs", MarkerSize=10, MarkerFaceColor="r")
    xlabel("X (mm)"); ylabel("Y (mm)")
    title("Pin A Path (Best Design)")
    legend("Path", "\theta_2 = 0°", "Goal", Location="best")
    grid on; axis equal; box on
    exportgraphics(gcf, fullfile(scriptDir, "ENME473_Q7_BestDesign.png"), Resolution=600);

    % plot best design every 15 deg from 0 to 120
    plotAnglesFigure(bx_best.allJoints, pA_best, 0:15:120, goalX, goalY, "Best");
    exportgraphics(gcf, fullfile(scriptDir, "ENME473_Q7_BestDesign_Angles.png"), Resolution=600);

    % summary table: box dimensions and area every 5 deg (best)
    printAngleTable(bx_best.allJoints, "Best", xlsxPath);

    % report second best (different parameter combination)
    if length(results) >= 2
        alt = results(2);
        fprintf("\n=== ALTERNATE DESIGN ===\n");
        fprintf("Ground pivot O4: (%.2f, %.2f) mm\n", alt.xO4, alt.yO4);
        fprintf("RA = %.2f mm,  R8 = %.2f mm\n", alt.RA, alt.R8);
        fprintf("Enclosing box: %.2f x %.2f mm,  Area = %.2f mm^2\n", ...
            alt.width, alt.height, alt.area);
        fprintf("Area reduction: %.2f%%\n", (1 - alt.area/origArea)*100);
        fprintf("Goal reached at theta2 = %d deg (distance = %.2f mm)\n", ...
            alt.goalTheta2, alt.goalDist);
    end
else
    fprintf("\nNo design found that reaches the goal with a smaller box.\n");
    fprintf("Reporting two alternate designs with their goal and box info.\n\n");

    % report two alternate designs that got close
    alt1_RA = RA_orig * 0.85;
    alt1_R8 = R8_orig;
    alt1_R1 = R1_orig;
    alt1_t1 = theta1_orig;
    [~, pA1, bx1] = solveDesign(R2, R23, R14, R3, R4, R36, R46, ...
        R5, R6, R7, alt1_R8, alt1_RA, alt1_R1, alt1_t1, ...
        alpha, beta, gamma, x0_base);
    dist1 = sqrt((pA1.Ax - goalX).^2 + (pA1.Ay - goalY).^2);
    [md1, mi1] = min(dist1);

    fprintf("=== ALTERNATE DESIGN 1 ===\n");
    fprintf("Ground pivot O4: (%.2f, %.2f) mm  (unchanged)\n", xO4_orig, yO4_orig);
    fprintf("RA = %.2f mm (reduced),  R8 = %.2f mm\n", alt1_RA, alt1_R8);
    fprintf("Enclosing box: %.2f x %.2f mm,  Area = %.2f mm^2\n", bx1.width, bx1.height, bx1.area);
    fprintf("Closest to goal at theta2 = %d deg: (%.2f, %.2f), dist = %.2f mm\n", ...
        mi1-1, pA1.Ax(mi1), pA1.Ay(mi1), md1);

    alt2_xO4 = xO4_orig + 20;
    alt2_yO4 = yO4_orig;
    alt2_R1 = sqrt(alt2_xO4^2 + alt2_yO4^2);
    alt2_t1 = atan2(alt2_yO4, alt2_xO4);
    [~, pA2, bx2] = solveDesign(R2, R23, R14, R3, R4, R36, R46, ...
        R5, R6, R7, R8_orig, RA_orig, alt2_R1, alt2_t1, ...
        alpha, beta, gamma, x0_base);
    dist2 = sqrt((pA2.Ax - goalX).^2 + (pA2.Ay - goalY).^2);
    [md2, mi2] = min(dist2);

    fprintf("\n=== ALTERNATE DESIGN 2 ===\n");
    fprintf("Ground pivot O4: (%.2f, %.2f) mm  (shifted +20 mm in x)\n", alt2_xO4, alt2_yO4);
    fprintf("RA = %.2f mm,  R8 = %.2f mm  (unchanged)\n", RA_orig, R8_orig);
    fprintf("Enclosing box: %.2f x %.2f mm,  Area = %.2f mm^2\n", bx2.width, bx2.height, bx2.area);
    fprintf("Closest to goal at theta2 = %d deg: (%.2f, %.2f), dist = %.2f mm\n", ...
        mi2-1, pA2.Ax(mi2), pA2.Ay(mi2), md2);
end

%% comparison table
fprintf("\n=== SUMMARY TABLE ===\n");
fprintf("%-20s %10s %10s %10s %10s %10s %12s\n", ...
    "Design", "xO4", "yO4", "RA", "R8", "Box Area", "Area Change");
fprintf("%-20s %10.1f %10.1f %10.1f %10.1f %10.0f %12s\n", ...
    "Original", xO4_orig, yO4_orig, RA_orig, R8_orig, origArea, "baseline");

if ~isempty(results)
    best = results(1);
    fprintf("%-20s %10.1f %10.1f %10.1f %10.1f %10.0f %11.1f%%\n", ...
        "Best", best.xO4, best.yO4, best.RA, best.R8, best.area, ...
        (1-best.area/origArea)*100);
    if length(results) >= 2
        alt = results(2);
        fprintf("%-20s %10.1f %10.1f %10.1f %10.1f %10.0f %11.1f%%\n", ...
            "Alternate", alt.xO4, alt.yO4, alt.RA, alt.R8, alt.area, ...
            (1-alt.area/origArea)*100);
    end
end

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [allAngles, pinA, boxInfo] = solveDesign(R2, R23, R14, R3, R4, ...
    R36, R46, R5, R6, R7, R8, RA, R1, theta1, alpha, beta, gamma, x0)
% solves the full kinematics for theta2 = 0:120 deg and returns
% posture angles, Pin A coordinates, and bounding box info at theta2 = 0.

    RAtot = R8 + RA;
    in = 0:120;
    theta2_vals = deg2rad(in);
    N = length(theta2_vals);

    % storage
    rawAngles = zeros(6, N);
    Ax = zeros(1, N);
    Ay = zeros(1, N);
    allJoints = zeros(11, 2, N);

    x = x0;
    maxIter = 200;
    tol = 1e-8;

    for k = 1:N
        theta2 = theta2_vals(k);

        for iter = 1:maxIter
            theta23 = x(1); theta14 = x(2); theta46 = x(3);
            theta36 = x(4); theta8  = x(5); theta7  = x(6);

            f = [R2*cos(theta2) + R23*cos(theta23) - R14*cos(theta14) - R1*cos(theta1);
                 R2*sin(theta2) + R23*sin(theta23) - R14*sin(theta14) - R1*sin(theta1);
                 R2*cos(theta2) + R4*cos(theta23+gamma) + R46*cos(theta46) - R36*cos(theta36) - R3*cos(theta14+alpha) - R1*cos(theta1);
                 R2*sin(theta2) + R4*sin(theta23+gamma) + R46*sin(theta46) - R36*sin(theta36) - R3*sin(theta14+alpha) - R1*sin(theta1);
                 R2*cos(theta2) + R4*cos(theta23+gamma) + R6*cos(theta46) - R8*cos(theta8) + R7*cos(theta7) - R5*cos(theta36+beta) - R3*cos(theta14+alpha) - R1*cos(theta1);
                 R2*sin(theta2) + R4*sin(theta23+gamma) + R6*sin(theta46) - R8*sin(theta8) + R7*sin(theta7) - R5*sin(theta36+beta) - R3*sin(theta14+alpha) - R1*sin(theta1)];

            J = [-R23*sin(theta23),  R14*sin(theta14),             0,                   0,              0,             0;
                  R23*cos(theta23), -R14*cos(theta14),             0,                   0,              0,             0;
                 -R4*sin(theta23+gamma),  R3*sin(theta14+alpha), -R46*sin(theta46),  R36*sin(theta36),  0,             0;
                  R4*cos(theta23+gamma), -R3*cos(theta14+alpha),  R46*cos(theta46), -R36*cos(theta36),  0,             0;
                 -R4*sin(theta23+gamma),  R3*sin(theta14+alpha),  -R6*sin(theta46),  R5*sin(theta36+beta), R8*sin(theta8), -R7*sin(theta7);
                  R4*cos(theta23+gamma), -R3*cos(theta14+alpha),   R6*cos(theta46), -R5*cos(theta36+beta),-R8*cos(theta8),  R7*cos(theta7)];

            dx = J \ f;
            x_new = x - dx;
            x_new = atan2(sin(x_new), cos(x_new));

            if norm(f, inf) < tol && norm(x_new - x, inf) < tol
                x = x_new;
                break;
            end
            x = x_new;
        end

        if iter == maxIter
            error("NR did not converge at theta2 = %g deg.", in(k));
        end

        rawAngles(:, k) = x;

        theta23 = x(1); theta14 = x(2); theta46 = x(3);
        theta36 = x(4); theta8 = x(5);  theta7 = x(6);

        % Pin A position
        Ax(k) = R2*cos(theta2) + R4*cos(theta23+gamma) + R6*cos(theta46) - RAtot*cos(theta8);
        Ay(k) = R2*sin(theta2) + R4*sin(theta23+gamma) + R6*sin(theta46) - RAtot*sin(theta8);

        % all joint positions at this theta2
        O2 = [0, 0];
        O4 = [R1*cos(theta1), R1*sin(theta1)];
        JointA  = [R2*cos(theta2), R2*sin(theta2)];
        Joint43 = JointA + [R23*cos(theta23), R23*sin(theta23)];
        JointB  = JointA + [R4*cos(theta23+gamma), R4*sin(theta23+gamma)];
        JointD  = O4 + [R3*cos(theta14+alpha), R3*sin(theta14+alpha)];
        JointC  = JointB + [R46*cos(theta46), R46*sin(theta46)];
        P1      = JointB + [R6*cos(theta46), R6*sin(theta46)];
        P2      = P1 - [R8*cos(theta8), R8*sin(theta8)];
        JointE7 = P2 + [R7*cos(theta7), R7*sin(theta7)];
        PinA_pos = [Ax(k), Ay(k)];
        allJoints(:,:,k) = [O2; O4; JointA; Joint43; JointB; JointD; JointC; P1; P2; JointE7; PinA_pos];
    end

    allAngles.raw = rawAngles;
    pinA.Ax = Ax;
    pinA.Ay = Ay;

    % bounding box at theta2 = 0
    allPts = allJoints(:,:,1);
    boxInfo = boxFromPts(allPts);
    boxInfo.joints = allPts;
    boxInfo.labels = ["O2"; "O4"; "A"; "43"; "B"; "D"; "C"; "P1"; "P2"; "E7"; "PinA"];
    boxInfo.allJoints = allJoints;
end

function bx = boxFromPts(pts)
% compute bounding box info from a set of (x,y) points
    bx.xmin   = min(pts(:,1));
    bx.xmax   = max(pts(:,1));
    bx.ymin   = min(pts(:,2));
    bx.ymax   = max(pts(:,2));
    bx.width  = bx.xmax - bx.xmin;
    bx.height = bx.ymax - bx.ymin;
    bx.area   = bx.width * bx.height;
end

function plotMechanism(joints, pinAx, pinAy, boxInfo)
% plots the mechanism skeleton and bounding box

    % joint coordinates
    O2 = joints(1,:);  O4 = joints(2,:);
    A  = joints(3,:);  J43 = joints(4,:);
    B  = joints(5,:);  D = joints(6,:);
    C  = joints(7,:);  P1 = joints(8,:);
    P2 = joints(9,:);  E7 = joints(10,:);

    hold on

    % draw links as lines
    % link 1 (ground)
    plot([O2(1) O4(1)], [O2(2) O4(2)], "k-", LineWidth=3)
    % link 2 (input)
    plot([O2(1) A(1)], [O2(2) A(2)], "b-", LineWidth=2)
    % link 3 (ternary: O4-J43-D)
    plot([O4(1) J43(1) D(1) O4(1)], [O4(2) J43(2) D(2) O4(2)], "r-", LineWidth=2)
    % link 4 (ternary: A-J43-B)
    plot([A(1) J43(1) B(1) A(1)], [A(2) J43(2) B(2) A(2)], "g-", LineWidth=2)
    % link 5 (ternary: D-C-E7)
    plot([D(1) C(1) E7(1) D(1)], [D(2) C(2) E7(2) D(2)], "m-", LineWidth=2)
    % link 6 (B-C-P1 collinear)
    plot([B(1) P1(1)], [B(2) P1(2)], Color=[0 0.6 0.6], LineWidth=2)
    % link 7
    plot([E7(1) P2(1)], [E7(2) P2(2)], Color=[0.6 0.3 0], LineWidth=2)
    % link 8 (P1-P2-PinA)
    plot([P1(1) P2(1) pinAx], [P1(2) P2(2) pinAy], Color=[0.5 0 0.5], LineWidth=2)

    % draw joints
    plot(joints(:,1), joints(:,2), "ko", MarkerSize=6, MarkerFaceColor="k")
    plot(pinAx, pinAy, "r^", MarkerSize=10, MarkerFaceColor="r")

    % draw bounding box
    rectangle(Position=[boxInfo.xmin, boxInfo.ymin, boxInfo.width, boxInfo.height], ...
        EdgeColor="r", LineStyle="--", LineWidth=1.5)

    % label Pin A
    text(pinAx + 15, pinAy, "Pin A", FontSize=9, FontWeight="bold")

    xlabel("X (mm)"); ylabel("Y (mm)")
    grid on; axis equal; box on
    hold off
end

function plotAnglesFigure(allJoints, pinA, anglesDeg, goalX, goalY, label)
% plot mechanism, bounding box, and pin A path at multiple theta2 angles
    n = length(anglesDeg);
    nCols = ceil(sqrt(n));
    nRows = ceil(n / nCols);
    figure("Position", [100, 100, 350*nCols, 300*nRows])
    for i = 1:n
        ang = anglesDeg(i);
        idx = ang + 1; % theta2 = 0:120 deg, 1-indexed
        subplot(nRows, nCols, i)
        jts = allJoints(:,:,idx);
        bx = boxFromPts(jts);
        plotMechanism(jts, pinA.Ax(idx), pinA.Ay(idx), bx);
        hold on
        plot(pinA.Ax(1:idx), pinA.Ay(1:idx), "b:", LineWidth=1.5)
        plot(goalX, goalY, "rs", MarkerSize=8, MarkerFaceColor="r")
        title(sprintf("\\theta_2 = %d°  (Box: %.0f × %.0f mm)", ang, bx.width, bx.height))
        hold off
    end
    sgtitle(label + " Design: Mechanism, Bounding Box, and Pin A Path")
end

function printAngleTable(allJoints, label, xlsxPath)
% print bounding box dimensions and area every 5 deg from 0 to 120
% and export the table to an Excel file (one sheet per design)
    angles = (0:5:120)';
    N = length(angles);
    width  = zeros(N, 1);
    height = zeros(N, 1);
    area   = zeros(N, 1);
    for i = 1:N
        idx = angles(i) + 1; % theta2 = 0:120 deg, 1-indexed
        bx = boxFromPts(allJoints(:,:,idx));
        width(i)  = bx.width;
        height(i) = bx.height;
        area(i)   = bx.area;
    end

    T = table(angles, width, height, area, ...
        VariableNames=["Theta2_deg", "Width_mm", "Height_mm", "Area_mm2"]);

    fprintf("\n=== %s DESIGN: BOX DIMENSIONS vs theta2 ===\n", upper(label));
    disp(T)

    writetable(T, xlsxPath, Sheet=label, WriteMode="overwritesheet");
    fprintf("Exported to: %s (sheet: %s)\n", xlsxPath, label);
end

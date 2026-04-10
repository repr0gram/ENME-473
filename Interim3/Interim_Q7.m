% ENME 473 Project - Part III: Design (Question 7)
% finds the "goal" Pin A position at theta2 = 120 deg,
% computes the enclosing box at theta2 = 0 deg,
% and systematically varies three design parameters:
%   (a) ground link angle theta1 (|R1| held fixed; right pin at origin)
%   (b) dP1A: distance from Pin A to the far link-8 pin  (originally 580.7)
%   (c) dP2A: distance from Pin A to the near link-8 pin (originally 336.9)
% Pin A is offset from the P1-P2 axis — its position is obtained by
% triangulating the rigid triangle formed by P1, P2, and Pin A on link 8.

clc; clear; close all;
scriptDir = fileparts(mfilename('fullpath'));

% single tabbed figure window to hold all four plots
mainFig = figure(Name="ENME 473 Q7 Results", Position=[100 100 1400 900]);
mainTabs = uitabgroup(mainFig);

% collect all mechanism-plot axes so they can share xlim/ylim at the end
mechAxes = gobjects(0);

%% original geometry (mm)
R1  = sqrt(237.2^2 + 70.6^2); % ground-link length is held FIXED across designs
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
R8  = 254.4;     % fixed: distance between the two pins on link 8

% link-8 triangle side lengths (design parameters b, c)
dP1A_orig = 580.7; % |Pin A - P1|, far link-8 pin
dP2A_orig = 336.9; % |Pin A - P2|, near link-8 pin

% which side of the P1->P2 axis Pin A sits on (+1 or -1)
% flip this if the plotted mechanism does not match the PDF figure
pinASide = -1;

theta1_orig = deg2rad(180 - atand(70.6/237.2));
alpha  = deg2rad(3.2);
beta   = deg2rad(3.5);
gamma  = deg2rad(1.8);

% ground pivot O4 original position (mm), derived from R1 and theta1
xO4_orig = R1*cos(theta1_orig);
yO4_orig = R1*sin(theta1_orig);

%% solve original design
fprintf("=== ORIGINAL DESIGN ===\n");

% initial guess for NR (converges to correct branch)
x0 = deg2rad([167.5; -4.0; 10.7; 165.4; -4.8; -171.3]);

% solve full range 0-120 deg
[allAngles, pinA, boxInfo] = solveDesign(R2, R23, R14, R3, R4, R36, R46, ...
    R5, R6, R7, R8, dP1A_orig, dP2A_orig, pinASide, R1, theta1_orig, ...
    alpha, beta, gamma, x0);

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
tabOrig = uitab(mainTabs, Title="Original");
tlOrig = tiledlayout(tabOrig, 1, 2);
nexttile(tlOrig)
plotMechanism(boxInfo.joints, pinA.Ax(1), pinA.Ay(1), boxInfo);
title("Mechanism at \theta_2 = 0° (Original)")
mechAxes(end+1) = gca;

nexttile(tlOrig)
plot(pinA.Ax, pinA.Ay, "b-", LineWidth=2)
hold on
plot(pinA.Ax(1), pinA.Ay(1), "go", MarkerSize=10, MarkerFaceColor="g")
plot(goalX, goalY, "rs", MarkerSize=10, MarkerFaceColor="r")
xlabel("X (mm)"); ylabel("Y (mm)")
title("Pin A Path (Original)")
legend("Path", "\theta_2 = 0°", "\theta_2 = 120° (Goal)", Location="best")
grid on; axis equal; box on
exportgraphics(tlOrig, fullfile(scriptDir, "ENME473_Q7_Original.png"), Resolution=600);

%% plot original mechanism and path every 15 deg from 0 to 120
tabOrigAng = uitab(mainTabs, Title="Original Angles");
tlOrigAng = plotAnglesFigure(boxInfo.allJoints, pinA, 0:15:120, goalX, goalY, "Original", tabOrigAng);
mechAxes = [mechAxes(:); findobj(tlOrigAng, Type="axes")];
exportgraphics(tlOrigAng, fullfile(scriptDir, "ENME473_Q7_Original_Angles.png"), Resolution=600);

%% summary table: box dimensions and area every 5 deg (original)
xlsxPath = fullfile(scriptDir, "ENME473_Q7_BoxDimensions.xlsx");
if isfile(xlsxPath)
    delete(xlsxPath) % start fresh each run
end
printAngleTable(boxInfo.allJoints, "Original", xlsxPath);

%% design parameter search
% three design parameters:
%   (a) theta1 — ground-link angle (|R1| fixed; O2 pinned at origin)
%   (b) dP1A  — distance from Pin A to the far link-8 pin (was 580.7)
%   (c) dP2A  — distance from Pin A to the near link-8 pin (was 336.9)
%
% for each candidate: solve 0..120 deg, check Pin A reaches near the goal,
% and keep those with a smaller bounding box at theta2 = 0.

goalTol = 15; % mm tolerance for reaching goal position

% search ranges
theta1_range = linspace(theta1_orig - deg2rad(15), theta1_orig + deg2rad(15), 21);
dP1A_range   = linspace(dP1A_orig*0.8, dP1A_orig*1.1, 16);
dP2A_range   = linspace(dP2A_orig*0.7, dP2A_orig*1.1, 17);

% converged angles from the original design — used as NR initial guess
x0_base = allAngles.raw(:, 1);

% storage for results
results = [];
nTotal = length(theta1_range) * length(dP1A_range) * length(dP2A_range);
fprintf("\n=== DESIGN SEARCH ===\n");
fprintf("Evaluating %d design candidates...\n", nTotal);

for iT = 1:length(theta1_range)
    theta1_c = theta1_range(iT);
    for iD1 = 1:length(dP1A_range)
        for iD2 = 1:length(dP2A_range)
            dP1A_c = dP1A_range(iD1);
            dP2A_c = dP2A_range(iD2);

            % triangle inequality must hold for a valid link-8 geometry
            if dP1A_c + dP2A_c <= R8 || abs(dP1A_c - dP2A_c) >= R8
                continue
            end

            try
                [~, pA, bx] = solveDesign(R2, R23, R14, R3, R4, R36, R46, ...
                    R5, R6, R7, R8, dP1A_c, dP2A_c, pinASide, ...
                    R1, theta1_c, alpha, beta, gamma, x0_base);

                % check if Pin A passes near goal at any theta2
                dist = sqrt((pA.Ax - goalX).^2 + (pA.Ay - goalY).^2);
                [minDist, idxMin] = min(dist);

                if minDist < goalTol && bx.area < origArea
                    r.theta1 = theta1_c;
                    r.dP1A   = dP1A_c;
                    r.dP2A   = dP2A_c;
                    r.xO4    = R1*cos(theta1_c);
                    r.yO4    = R1*sin(theta1_c);
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
    fprintf("Ground pivot O4: (%.2f, %.2f) mm  [theta1 = %.2f deg]\n", ...
        best.xO4, best.yO4, rad2deg(best.theta1));
    fprintf("dP1A = %.2f mm  (was %.2f)\n", best.dP1A, dP1A_orig);
    fprintf("dP2A = %.2f mm  (was %.2f)\n", best.dP2A, dP2A_orig);
    fprintf("Enclosing box: %.2f x %.2f mm,  Area = %.2f mm^2\n", ...
        best.width, best.height, best.area);
    fprintf("Area reduction: %.2f%%\n", (1 - best.area/origArea)*100);
    fprintf("Goal reached at theta2 = %d deg (distance = %.2f mm)\n", ...
        best.goalTheta2, best.goalDist);
    fprintf("Pin A at goal: (%.2f, %.2f) mm\n", best.goalAx, best.goalAy);

    % plot best design
    [~, pA_best, bx_best] = solveDesign(R2, R23, R14, R3, R4, R36, R46, ...
        R5, R6, R7, R8, best.dP1A, best.dP2A, pinASide, ...
        R1, best.theta1, alpha, beta, gamma, x0_base);

    tabBest = uitab(mainTabs, Title="Best");
    tlBest = tiledlayout(tabBest, 1, 2);
    nexttile(tlBest)
    plotMechanism(bx_best.joints, pA_best.Ax(1), pA_best.Ay(1), bx_best);
    title("Mechanism at \theta_2 = 0° (Best Design)")
    mechAxes(end+1) = gca;

    nexttile(tlBest)
    plot(pA_best.Ax, pA_best.Ay, "b-", LineWidth=2)
    hold on
    plot(pA_best.Ax(1), pA_best.Ay(1), "go", MarkerSize=10, MarkerFaceColor="g")
    plot(goalX, goalY, "rs", MarkerSize=10, MarkerFaceColor="r")
    xlabel("X (mm)"); ylabel("Y (mm)")
    title("Pin A Path (Best Design)")
    legend("Path", "\theta_2 = 0°", "Goal", Location="best")
    grid on; axis equal; box on
    exportgraphics(tlBest, fullfile(scriptDir, "ENME473_Q7_BestDesign.png"), Resolution=600);

    % plot best design every 15 deg from 0 to 120
    tabBestAng = uitab(mainTabs, Title="Best Angles");
    tlBestAng = plotAnglesFigure(bx_best.allJoints, pA_best, 0:15:120, goalX, goalY, "Best", tabBestAng);
    mechAxes = [mechAxes(:); findobj(tlBestAng, Type="axes")];
    exportgraphics(tlBestAng, fullfile(scriptDir, "ENME473_Q7_BestDesign_Angles.png"), Resolution=600);

    % unify xlim/ylim across all mechanism plots so original vs best
    % are directly comparable. pad by 5% of each axis range.
    allPts = cat(1, boxInfo.allJoints, bx_best.allJoints);
    xMin = min(allPts(:,1,:), [], "all");
    xMax = max(allPts(:,1,:), [], "all");
    yMin = min(allPts(:,2,:), [], "all");
    yMax = max(allPts(:,2,:), [], "all");
    padX = 0.05 * (xMax - xMin);
    padY = 0.05 * (yMax - yMin);
    sharedXLim = [xMin - padX, xMax + padX];
    sharedYLim = [yMin - padY, yMax + padY];
    for i = 1:numel(mechAxes)
        xlim(mechAxes(i), sharedXLim)
        ylim(mechAxes(i), sharedYLim)
    end

    % re-export figures now that limits have been unified
    exportgraphics(tlOrig, fullfile(scriptDir, "ENME473_Q7_Original.png"), Resolution=600);
    exportgraphics(tlOrigAng, fullfile(scriptDir, "ENME473_Q7_Original_Angles.png"), Resolution=600);
    exportgraphics(tlBest, fullfile(scriptDir, "ENME473_Q7_BestDesign.png"), Resolution=600);
    exportgraphics(tlBestAng, fullfile(scriptDir, "ENME473_Q7_BestDesign_Angles.png"), Resolution=600);

    % summary table: box dimensions and area every 5 deg (best)
    printAngleTable(bx_best.allJoints, "Best", xlsxPath);

    % report second best (different parameter combination)
    if length(results) >= 2
        alt = results(2);
        fprintf("\n=== ALTERNATE DESIGN ===\n");
        fprintf("Ground pivot O4: (%.2f, %.2f) mm  [theta1 = %.2f deg]\n", ...
            alt.xO4, alt.yO4, rad2deg(alt.theta1));
        fprintf("dP1A = %.2f mm,  dP2A = %.2f mm\n", alt.dP1A, alt.dP2A);
        fprintf("Enclosing box: %.2f x %.2f mm,  Area = %.2f mm^2\n", ...
            alt.width, alt.height, alt.area);
        fprintf("Area reduction: %.2f%%\n", (1 - alt.area/origArea)*100);
        fprintf("Goal reached at theta2 = %d deg (distance = %.2f mm)\n", ...
            alt.goalTheta2, alt.goalDist);
    end
else
    fprintf("\nNo design found that reaches the goal with a smaller box.\n");
    fprintf("Reporting two alternate designs with their goal and box info.\n\n");

    % alternate 1: shrink dP1A by 15%
    alt1_dP1A = dP1A_orig * 0.85;
    alt1_dP2A = dP2A_orig;
    alt1_t1   = theta1_orig;
    [~, pA1, bx1] = solveDesign(R2, R23, R14, R3, R4, R36, R46, ...
        R5, R6, R7, R8, alt1_dP1A, alt1_dP2A, pinASide, ...
        R1, alt1_t1, alpha, beta, gamma, x0_base);
    dist1 = sqrt((pA1.Ax - goalX).^2 + (pA1.Ay - goalY).^2);
    [md1, mi1] = min(dist1);

    fprintf("=== ALTERNATE DESIGN 1 ===\n");
    fprintf("Ground pivot O4: (%.2f, %.2f) mm  (unchanged)\n", xO4_orig, yO4_orig);
    fprintf("dP1A = %.2f mm (reduced),  dP2A = %.2f mm\n", alt1_dP1A, alt1_dP2A);
    fprintf("Enclosing box: %.2f x %.2f mm,  Area = %.2f mm^2\n", bx1.width, bx1.height, bx1.area);
    fprintf("Closest to goal at theta2 = %d deg: (%.2f, %.2f), dist = %.2f mm\n", ...
        mi1-1, pA1.Ax(mi1), pA1.Ay(mi1), md1);

    % alternate 2: rotate the ground link by +5 deg
    alt2_t1   = theta1_orig + deg2rad(5);
    alt2_dP1A = dP1A_orig;
    alt2_dP2A = dP2A_orig;
    [~, pA2, bx2] = solveDesign(R2, R23, R14, R3, R4, R36, R46, ...
        R5, R6, R7, R8, alt2_dP1A, alt2_dP2A, pinASide, ...
        R1, alt2_t1, alpha, beta, gamma, x0_base);
    dist2 = sqrt((pA2.Ax - goalX).^2 + (pA2.Ay - goalY).^2);
    [md2, mi2] = min(dist2);

    fprintf("\n=== ALTERNATE DESIGN 2 ===\n");
    fprintf("Ground pivot O4: (%.2f, %.2f) mm  [theta1 = %.2f deg, rotated +5 deg]\n", ...
        R1*cos(alt2_t1), R1*sin(alt2_t1), rad2deg(alt2_t1));
    fprintf("dP1A = %.2f mm,  dP2A = %.2f mm  (unchanged)\n", alt2_dP1A, alt2_dP2A);
    fprintf("Enclosing box: %.2f x %.2f mm,  Area = %.2f mm^2\n", bx2.width, bx2.height, bx2.area);
    fprintf("Closest to goal at theta2 = %d deg: (%.2f, %.2f), dist = %.2f mm\n", ...
        mi2-1, pA2.Ax(mi2), pA2.Ay(mi2), md2);
end

%% comparison table
fprintf("\n=== SUMMARY TABLE ===\n");
fprintf("%-12s %10s %10s %10s %10s %10s %12s\n", ...
    "Design", "xO4", "yO4", "dP1A", "dP2A", "Box Area", "Area Change");
fprintf("%-12s %10.1f %10.1f %10.1f %10.1f %10.0f %12s\n", ...
    "Original", xO4_orig, yO4_orig, dP1A_orig, dP2A_orig, origArea, "baseline");

if ~isempty(results)
    best = results(1);
    fprintf("%-12s %10.1f %10.1f %10.1f %10.1f %10.0f %11.1f%%\n", ...
        "Best", best.xO4, best.yO4, best.dP1A, best.dP2A, best.area, ...
        (1-best.area/origArea)*100);
    if length(results) >= 2
        alt = results(2);
        fprintf("%-12s %10.1f %10.1f %10.1f %10.1f %10.0f %11.1f%%\n", ...
            "Alternate", alt.xO4, alt.yO4, alt.dP1A, alt.dP2A, alt.area, ...
            (1-alt.area/origArea)*100);
    end
end

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [allAngles, pinA, boxInfo] = solveDesign(R2, R23, R14, R3, R4, ...
    R36, R46, R5, R6, R7, R8, dP1A, dP2A, pinASide, R1, theta1, ...
    alpha, beta, gamma, x0)
% solves the full kinematics for theta2 = 0:120 deg and returns
% posture angles, Pin A coordinates, and bounding box info at theta2 = 0.
% Pin A is NOT collinear with the two link-8 pins — its position is
% triangulated from the rigid triangle P1-P2-PinA with sides R8, dP1A, dP2A.

    % local-frame (x along P1->P2) coordinates of Pin A
    pinA_xL = (dP1A^2 - dP2A^2 + R8^2) / (2*R8);
    pinA_yL = pinASide * sqrt(max(0, dP1A^2 - pinA_xL^2));

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
        theta8 = x(5);  theta7 = x(6);

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

        % triangulate Pin A from the rigid P1-P2-PinA triangle
        % unit vector along P1->P2 and its CCW perpendicular
        u = (P2 - P1) / R8;
        v = [-u(2), u(1)];
        PinA_pos = P1 + pinA_xL*u + pinA_yL*v;
        Ax(k) = PinA_pos(1);
        Ay(k) = PinA_pos(2);

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
    % link 8 (P1-P2-PinA triangle)
    plot([P1(1) P2(1) pinAx P1(1)], [P1(2) P2(2) pinAy P1(2)], Color=[0.5 0 0.5], LineWidth=2)

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

function tl = plotAnglesFigure(allJoints, pinA, anglesDeg, goalX, goalY, label, parent)
% plot mechanism, bounding box, and pin A path at multiple theta2 angles
    n = length(anglesDeg);
    nCols = ceil(sqrt(n));
    nRows = ceil(n / nCols);
    tl = tiledlayout(parent, nRows, nCols);
    for i = 1:n
        ang = anglesDeg(i);
        idx = ang + 1; % theta2 = 0:120 deg, 1-indexed
        nexttile(tl)
        jts = allJoints(:,:,idx);
        bx = boxFromPts(jts);
        plotMechanism(jts, pinA.Ax(idx), pinA.Ay(idx), bx);
        hold on
        plot(pinA.Ax(1:idx), pinA.Ay(1:idx), "b:", LineWidth=1.5)
        plot(goalX, goalY, "rs", MarkerSize=8, MarkerFaceColor="r")
        title(sprintf("\\theta_2 = %d°  (Box: %.0f × %.0f mm)", ang, bx.width, bx.height))
        hold off
    end
    title(tl, label + " Design: Mechanism, Bounding Box, and Pin A Path")
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

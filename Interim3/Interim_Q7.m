% ENME 473 Project - Part III: Design (Question 7)
% solves the Q7 design problem "backwards": instead of sweeping candidate
% (theta1, dP1A, dP2A) triples and checking whether Pin A happens to pass
% near the goal, this script fixes Pin A = goal EXACTLY and computes the
% required dP1A, dP2A analytically:
%   1. for each candidate theta1 and each target input angle theta2_target,
%      solve the mechanism to get P1 and P2 (these depend only on theta1,
%      not on the link-8 triangle sides)
%   2. compute dP1A = |P1 - goal| and dP2A = |P2 - goal|
%   3. determine pinASide from which side of the P1->P2 axis the goal lies
%   4. re-place Pin A on the joint tensor and evaluate the bounding box at
%      theta2 = 0
% every accepted design therefore hits the goal exactly (to machine
% precision), and the entire search budget is spent on minimizing the
% bounding box rather than hunting for goal-reach.

clc; clear; close all;
scriptDir = fileparts(mfilename('fullpath'));

% single tabbed figure window to hold all plots
mainFig = figure(Name="ENME 473 Q7 Results", Position=[100 100 1400 900]);
mainTabs = uitabgroup(mainFig);

% single-frame mechanism plots (Original + Best tabs) and multi-frame
% mechanism grids (Original Angles + Best Angles tabs) are tracked
% separately so they can be given different shared limits: the single
% frames zoom tight to theta2 = 0, the angle grids span the full sweep.
mechAxesSingle = gobjects(0);
mechAxesAngles = gobjects(0);

%% original geometry (mm)
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

dP1A_orig = 580.7;
dP2A_orig = 336.9;
pinASide_orig = -1;

theta1_orig = deg2rad(180 - atand(70.6/237.2));
alpha  = deg2rad(3.2);
beta   = deg2rad(3.5);
gamma  = deg2rad(1.8);

xO4_orig = R1*cos(theta1_orig);
yO4_orig = R1*sin(theta1_orig);

%% solve original design (baseline)
fprintf("=== ORIGINAL DESIGN ===\n");

x0 = deg2rad([167.5; -4.0; 10.7; 165.4; -4.8; -171.3]);

[allAngles, pinA, boxInfo] = solveDesign(R2, R23, R14, R3, R4, R36, R46, ...
    R5, R6, R7, R8, dP1A_orig, dP2A_orig, pinASide_orig, R1, theta1_orig, ...
    alpha, beta, gamma, x0);

goalX = pinA.Ax(end);
goalY = pinA.Ay(end);
goal  = [goalX, goalY];
fprintf("Pin A at theta2 = 120 deg (GOAL): (%.2f, %.2f) mm\n", goalX, goalY);

fprintf("Enclosing box at theta2 = 0 deg:\n");
fprintf("  Width  = %.2f mm\n", boxInfo.width);
fprintf("  Height = %.2f mm\n", boxInfo.height);
fprintf("  Area   = %.2f mm^2\n", boxInfo.area);

origArea = boxInfo.area;

%% plot original
tabOrig = uitab(mainTabs, Title="Original");
axOrig = axes(tabOrig);
plotMechanism(boxInfo.joints, pinA.Ax(1), pinA.Ay(1), boxInfo);
title(axOrig, "Mechanism at \theta_2 = 0° (Original)")
mechAxesSingle(end+1) = axOrig;
exportgraphics(axOrig, fullfile(scriptDir, "ENME473_Q7_Original.png"), Resolution=600);

tabOrigAng = uitab(mainTabs, Title="Original Angles");
tlOrigAng = plotAnglesFigure(boxInfo.allJoints, pinA, 0:15:120, goalX, goalY, "Original", tabOrigAng);
mechAxesAngles = [mechAxesAngles(:); findobj(tlOrigAng, Type="axes")];
exportgraphics(tlOrigAng, fullfile(scriptDir, "ENME473_Q7_Original_Angles.png"), Resolution=600);

%% original box table (every 1 deg)
xlsxPath = fullfile(scriptDir, "ENME473_Q7_BoxDimensions.xlsx");
if isfile(xlsxPath)
    delete(xlsxPath)
end
printAngleTable(boxInfo.allJoints, "Original", xlsxPath);

%% design search
% only theta1 and the target input angle are swept. dP1A and dP2A are
% computed analytically from the goal, so every surviving candidate has
% Pin A = goal at theta2 = theta2_target, exactly.

fprintf("\n=== DESIGN SEARCH ===\n");

theta1_range = linspace(theta1_orig - deg2rad(20), theta1_orig + deg2rad(20), 81);
theta2_target_range = 90:120;

% keep the computed link-8 triangle sides within a sensible range so we
% do not report designs that would require a dramatically different link 8
dP1A_bounds = [dP1A_orig*0.5, dP1A_orig*1.5];
dP2A_bounds = [dP2A_orig*0.5, dP2A_orig*1.5];

x0_base = allAngles.raw(:, 1);

nT  = length(theta1_range);
nTA = length(theta2_target_range);
fprintf("Evaluating %d theta1 values x %d target angles = %d candidates...\n", ...
    nT, nTA, nT*nTA);

results = [];

for iT = 1:nT
    theta1_c = theta1_range(iT);

    % solve the mechanism once for this theta1. dP1A/dP2A passed in here
    % are placeholders — the kinematics do not depend on them, and we will
    % re-place Pin A below using the values required to hit the goal.
    try
        [~, ~, bxTmp] = solveDesign(R2, R23, R14, R3, R4, R36, R46, ...
            R5, R6, R7, R8, dP1A_orig, dP2A_orig, pinASide_orig, ...
            R1, theta1_c, alpha, beta, gamma, x0_base);
    catch
        continue
    end

    jointsRaw = bxTmp.allJoints; % 11 x 2 x 121

    for iA = 1:nTA
        theta2_target = theta2_target_range(iA);
        idx = theta2_target + 1; % theta2 = 0..120 deg, 1-indexed

        P1 = squeeze(jointsRaw(8, :, idx));
        P2 = squeeze(jointsRaw(9, :, idx));

        dP1A_c = norm(goal - P1);
        dP2A_c = norm(goal - P2);

        % physical plausibility bounds on the link-8 triangle sides
        if dP1A_c < dP1A_bounds(1) || dP1A_c > dP1A_bounds(2)
            continue
        end
        if dP2A_c < dP2A_bounds(1) || dP2A_c > dP2A_bounds(2)
            continue
        end

        % triangle inequality on (R8, dP1A, dP2A)
        if dP1A_c + dP2A_c <= R8 || abs(dP1A_c - dP2A_c) >= R8
            continue
        end

        % which side of the P1->P2 axis the goal lies on
        u  = (P2 - P1) / R8;
        v  = [-u(2), u(1)];
        yL = dot(goal - P1, v);
        if abs(yL) < 1e-9
            continue % goal collinear with P1-P2, degenerate triangle
        end
        pinASide_c = sign(yL);

        % re-place Pin A at every theta2 using the required triangle sides
        newJoints = replacePinA(jointsRaw, dP1A_c, dP2A_c, pinASide_c, R8);

        bxNew = boxFromPts(newJoints(:,:,1));

        if bxNew.area >= origArea
            continue
        end

        r.theta1        = theta1_c;
        r.theta2_target = theta2_target;
        r.dP1A          = dP1A_c;
        r.dP2A          = dP2A_c;
        r.pinASide      = pinASide_c;
        r.xO4           = R1*cos(theta1_c);
        r.yO4           = R1*sin(theta1_c);
        r.area          = bxNew.area;
        r.width         = bxNew.width;
        r.height        = bxNew.height;
        results = [results; r]; %#ok<AGROW>
    end
end

fprintf("Search complete. Found %d designs with smaller box and EXACT goal match.\n", ...
    length(results));

%% report results
if ~isempty(results)
    areas = [results.area];
    [~, sortIdx] = sort(areas);
    results = results(sortIdx);

    best = results(1);
    fprintf("\n=== BEST DESIGN ===\n");
    fprintf("Ground pivot O4: (%.2f, %.2f) mm  [theta1 = %.2f deg]\n", ...
        best.xO4, best.yO4, rad2deg(best.theta1));
    fprintf("dP1A = %.2f mm  (was %.2f)\n", best.dP1A, dP1A_orig);
    fprintf("dP2A = %.2f mm  (was %.2f)\n", best.dP2A, dP2A_orig);
    fprintf("Goal reached EXACTLY at theta2 = %d deg\n", best.theta2_target);
    fprintf("Enclosing box: %.2f x %.2f mm,  Area = %.2f mm^2\n", ...
        best.width, best.height, best.area);
    fprintf("Area reduction: %.2f%%\n", (1 - best.area/origArea)*100);

    % re-solve and plot
    [~, pA_best, bx_best] = solveDesign(R2, R23, R14, R3, R4, R36, R46, ...
        R5, R6, R7, R8, best.dP1A, best.dP2A, best.pinASide, ...
        R1, best.theta1, alpha, beta, gamma, x0_base);

    % sanity check: Pin A at theta2_target must equal the goal
    gotX = pA_best.Ax(best.theta2_target + 1);
    gotY = pA_best.Ay(best.theta2_target + 1);
    err  = norm([gotX - goalX, gotY - goalY]);
    fprintf("Pin A at theta2 = %d deg: (%.4f, %.4f)  (goal: %.4f, %.4f)\n", ...
        best.theta2_target, gotX, gotY, goalX, goalY);
    fprintf("Residual error: %.2e mm\n", err);

    tabBest = uitab(mainTabs, Title="Best");
    axBest = axes(tabBest);
    plotMechanism(bx_best.joints, pA_best.Ax(1), pA_best.Ay(1), bx_best);
    title(axBest, "Mechanism at \theta_2 = 0° (Best Design)")
    mechAxesSingle(end+1) = axBest;
    exportgraphics(axBest, fullfile(scriptDir, "ENME473_Q7_BestDesign.png"), Resolution=600);

    tabBestAng = uitab(mainTabs, Title="Best Angles");
    tlBestAng = plotAnglesFigure(bx_best.allJoints, pA_best, 0:15:120, goalX, goalY, "Best", tabBestAng);
    mechAxesAngles = [mechAxesAngles(:); findobj(tlBestAng, Type="axes")];
    exportgraphics(tlBestAng, fullfile(scriptDir, "ENME473_Q7_BestDesign_Angles.png"), Resolution=600);

    % Pin A paths overlay: both trajectories on a single axis for comparison
    tabPaths = uitab(mainTabs, Title="Pin A Paths");
    axPaths = axes(tabPaths);
    hold(axPaths, "on")
    plot(axPaths, pinA.Ax,  pinA.Ay,  "-", Color="b",           LineWidth=2, DisplayName="Original path")
    plot(axPaths, pA_best.Ax, pA_best.Ay, "-", Color=[0 0.6 0], LineWidth=2, DisplayName="Best path")
    plot(axPaths, pinA.Ax(1),  pinA.Ay(1),  "o", Color="b",     MarkerSize=10, MarkerFaceColor="b", DisplayName="Original \theta_2 = 0°")
    plot(axPaths, pA_best.Ax(1), pA_best.Ay(1), "o", Color=[0 0.6 0], MarkerSize=10, MarkerFaceColor=[0 0.6 0], DisplayName="Best \theta_2 = 0°")
    plot(axPaths, goalX, goalY, "rs", MarkerSize=12, MarkerFaceColor="r", DisplayName="Goal")
    hold(axPaths, "off")
    xlabel(axPaths, "X (mm)"); ylabel(axPaths, "Y (mm)")
    title(axPaths, "Pin A Paths: Original vs Best")
    legend(axPaths, Location="best")
    grid(axPaths, "on"); axis(axPaths, "equal"); box(axPaths, "on")
    exportgraphics(axPaths, fullfile(scriptDir, "ENME473_Q7_PinAPaths.png"), Resolution=600);

    % single-frame tabs (Original, Best): zoom tight to theta2 = 0 joints of
    % both designs with a small padding, and share the same limits so the two
    % mechanisms render at the same scale for direct comparison.
    singlePts = [boxInfo.joints; bx_best.joints];
    sxMin = min(singlePts(:,1));
    sxMax = max(singlePts(:,1));
    syMin = min(singlePts(:,2));
    syMax = max(singlePts(:,2));
    sPadX = 0.05 * (sxMax - sxMin);
    sPadY = 0.05 * (syMax - syMin);
    singleXLim = [sxMin - sPadX, sxMax + sPadX];
    singleYLim = [syMin - sPadY, syMax + sPadY];
    for i = 1:numel(mechAxesSingle)
        xlim(mechAxesSingle(i), singleXLim)
        ylim(mechAxesSingle(i), singleYLim)
    end

    % angle-grid tabs (Original Angles, Best Angles): wide shared limits that
    % accommodate the full theta2 = 0..120 sweep of both designs.
    allPts = cat(1, boxInfo.allJoints, bx_best.allJoints);
    axMin = min(allPts(:,1,:), [], "all");
    axMax = max(allPts(:,1,:), [], "all");
    ayMin = min(allPts(:,2,:), [], "all");
    ayMax = max(allPts(:,2,:), [], "all");
    aPadX = 0.05 * (axMax - axMin);
    aPadY = 0.05 * (ayMax - ayMin);
    anglesXLim = [axMin - aPadX, axMax + aPadX];
    anglesYLim = [ayMin - aPadY, ayMax + aPadY];
    for i = 1:numel(mechAxesAngles)
        xlim(mechAxesAngles(i), anglesXLim)
        ylim(mechAxesAngles(i), anglesYLim)
    end

    exportgraphics(axOrig, fullfile(scriptDir, "ENME473_Q7_Original.png"), Resolution=600);
    exportgraphics(tlOrigAng, fullfile(scriptDir, "ENME473_Q7_Original_Angles.png"), Resolution=600);
    exportgraphics(axBest, fullfile(scriptDir, "ENME473_Q7_BestDesign.png"), Resolution=600);
    exportgraphics(tlBestAng, fullfile(scriptDir, "ENME473_Q7_BestDesign_Angles.png"), Resolution=600);

    printAngleTable(bx_best.allJoints, "Best", xlsxPath);

    if length(results) >= 2
        alt = results(2);
        fprintf("\n=== ALTERNATE DESIGN ===\n");
        fprintf("Ground pivot O4: (%.2f, %.2f) mm  [theta1 = %.2f deg]\n", ...
            alt.xO4, alt.yO4, rad2deg(alt.theta1));
        fprintf("dP1A = %.2f mm,  dP2A = %.2f mm\n", alt.dP1A, alt.dP2A);
        fprintf("Goal reached EXACTLY at theta2 = %d deg\n", alt.theta2_target);
        fprintf("Enclosing box: %.2f x %.2f mm,  Area = %.2f mm^2\n", ...
            alt.width, alt.height, alt.area);
        fprintf("Area reduction: %.2f%%\n", (1 - alt.area/origArea)*100);
    end
else
    fprintf("\nNo inverse design found with a smaller box at theta2 = 0.\n");
end

%% summary table
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

%  ========================================================================
%%  LOCAL FUNCTIONS
%  ========================================================================

function newJoints = replacePinA(joints, dP1A, dP2A, pinASide, R8)
% re-place Pin A (joint index 11) at every theta2 frame of the joint tensor
% using the rigid P1-P2-PinA triangle with the given side lengths.
    pinA_xL = (dP1A^2 - dP2A^2 + R8^2) / (2*R8);
    pinA_yL = pinASide * sqrt(max(0, dP1A^2 - pinA_xL^2));

    N = size(joints, 3);
    newJoints = joints;
    for k = 1:N
        P1 = joints(8, :, k);
        P2 = joints(9, :, k);
        u = (P2 - P1) / R8;
        v = [-u(2), u(1)];
        newJoints(11, :, k) = P1 + pinA_xL*u + pinA_yL*v;
    end
end

function [allAngles, pinA, boxInfo] = solveDesign(R2, R23, R14, R3, R4, ...
    R36, R46, R5, R6, R7, R8, dP1A, dP2A, pinASide, R1, theta1, ...
    alpha, beta, gamma, x0)
% solves the full kinematics for theta2 = 0:120 deg and returns posture
% angles, Pin A coordinates, and bounding box info at theta2 = 0.

    pinA_xL = (dP1A^2 - dP2A^2 + R8^2) / (2*R8);
    pinA_yL = pinASide * sqrt(max(0, dP1A^2 - pinA_xL^2));

    in = 0:120;
    theta2_vals = deg2rad(in);
    N = length(theta2_vals);

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

    allPts = allJoints(:,:,1);
    boxInfo = boxFromPts(allPts);
    boxInfo.joints = allPts;
    boxInfo.labels = ["O2"; "O4"; "A"; "43"; "B"; "D"; "C"; "P1"; "P2"; "E7"; "PinA"];
    boxInfo.allJoints = allJoints;
end

function bx = boxFromPts(pts)
% axis-aligned bounding box info from a set of (x,y) points
    bx.xmin   = min(pts(:,1));
    bx.xmax   = max(pts(:,1));
    bx.ymin   = min(pts(:,2));
    bx.ymax   = max(pts(:,2));
    bx.width  = bx.xmax - bx.xmin;
    bx.height = bx.ymax - bx.ymin;
    bx.area   = bx.width * bx.height;
end

function plotMechanism(joints, pinAx, pinAy, boxInfo)
% draws the mechanism skeleton and bounding box

    O2 = joints(1,:);  O4 = joints(2,:);
    A  = joints(3,:);  J43 = joints(4,:);
    B  = joints(5,:);  D = joints(6,:);
    C  = joints(7,:);  P1 = joints(8,:);
    P2 = joints(9,:);  E7 = joints(10,:);

    hold on

    plot([O2(1) O4(1)], [O2(2) O4(2)], "k-", LineWidth=3)
    plot([O2(1) A(1)], [O2(2) A(2)], "b-", LineWidth=2)
    plot([O4(1) J43(1) D(1) O4(1)], [O4(2) J43(2) D(2) O4(2)], "r-", LineWidth=2)
    plot([A(1) J43(1) B(1) A(1)], [A(2) J43(2) B(2) A(2)], "g-", LineWidth=2)
    plot([D(1) C(1) E7(1) D(1)], [D(2) C(2) E7(2) D(2)], "m-", LineWidth=2)
    plot([B(1) P1(1)], [B(2) P1(2)], Color=[0 0.6 0.6], LineWidth=2)
    plot([E7(1) P2(1)], [E7(2) P2(2)], Color=[0.6 0.3 0], LineWidth=2)
    plot([P1(1) P2(1) pinAx P1(1)], [P1(2) P2(2) pinAy P1(2)], Color=[0.5 0 0.5], LineWidth=2)

    plot(joints(:,1), joints(:,2), "ko", MarkerSize=6, MarkerFaceColor="k")
    plot(pinAx, pinAy, "r^", MarkerSize=10, MarkerFaceColor="r")

    rectangle(Position=[boxInfo.xmin, boxInfo.ymin, boxInfo.width, boxInfo.height], ...
        EdgeColor="r", LineStyle="--", LineWidth=1.5)

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
        idx = ang + 1;
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
% print bounding box dimensions and area every 1 deg from 0 to 120
    angles = (0:1:120)';
    N = length(angles);
    width  = zeros(N, 1);
    height = zeros(N, 1);
    area   = zeros(N, 1);
    for i = 1:N
        idx = angles(i) + 1;
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

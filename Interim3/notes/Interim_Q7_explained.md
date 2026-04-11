# Interim_Q7.m — Walkthrough

This note explains `Interim3/Interim_Q7.m` end-to-end so someone unfamiliar with the
script can reproduce the workflow: solve the original 8-link mechanism, search a
3-parameter design space for a smaller bounding box, and produce plots + an Excel
report.

---

## 1. What the script is doing (high level)

The mechanism is an 8-link planar linkage driven by an input crank `theta2`. We care
about the trajectory of **Pin A** — a point fixed on the rigid triangle of link 8.
The goal is:

1. Solve the kinematics for `theta2 = 0..120°` with the original geometry.
2. Record the **goal**: where Pin A sits at `theta2 = 120°`.
3. Compute the **enclosing (axis-aligned) bounding box** of the whole mechanism at
   `theta2 = 0°`.
4. Sweep three design parameters and find a configuration whose mechanism has a
   smaller box at `theta2 = 0°` *and* still passes through (within a tolerance of)
   the original goal at some `theta2`.
5. Plot original vs. best, write a per-5° box-dimension table to Excel.

The three design knobs are:
- `theta1` — angle of the ground link (`|R1|` is held fixed; `O2` stays at the origin)
- `dP1A` — distance from Pin A to the *far* link-8 pin (P1)
- `dP2A` — distance from Pin A to the *near* link-8 pin (P2)

`R8` (the distance between the two link-8 pins themselves) is fixed.

---

## 2. Mechanism diagram

```
                        link 5 (D-C-E7 ternary)
                              D------C------E7
                             /              \
                            /                \  link 7
                           /                  \
                          /                    P2-----------+
                         /                     |            |
                        /                  link 8 (P1-P2-PinA rigid triangle)
              link 3   /                       |            |
              (O4-J43-D ternary)               P1           PinA  <-- tracked point
                      /                        |
                     /                         |
                    O4                         |
                    |  ground link (R1)        | link 6 (B-C-P1 collinear)
                    |  fixed length            |
                    |  angle = theta1          |
                    |                          B
                    |                         /|
                    |                        / |
                    |          link 4       /  |
                    |          (A-J43-B    /   |
                    |           ternary)  J43  |
                    |                    /     |
                    O2 ---------------- A      |
                       link 2 (input)    \____/
                       angle = theta2
```

Joint name reference (matches the `joints` array indices in the code):

| idx | name  | description                                 |
|-----|-------|---------------------------------------------|
| 1   | O2    | input pivot, fixed at origin                |
| 2   | O4    | ground pivot, at `R1·(cos θ1, sin θ1)`      |
| 3   | A     | end of link 2 (input crank tip)             |
| 4   | J43   | shared corner of link-3 and link-4 ternaries|
| 5   | B     | far corner of link-4 ternary                |
| 6   | D     | far corner of link-3 ternary                |
| 7   | C     | shared corner of link-5 / link-6 path       |
| 8   | P1    | far link-8 pin                              |
| 9   | P2    | near link-8 pin                             |
| 10  | E7    | end of link 7                               |
| 11  | PinA  | tracked point on link 8 (NOT collinear)     |

Pin A is **off** the P1–P2 line; its position is solved by triangulating the rigid
triangle whose sides are `R8`, `dP1A`, and `dP2A`. The sign `pinASide = ±1` picks
which side of the P1→P2 axis Pin A sits on.

---

## 3. File / code structure

```
Interim_Q7.m
├── (script body)
│   ├── §1  fixed link lengths + nominal angles (alpha, beta, gamma, theta1)
│   ├── §2  solve original design (1 call to solveDesign)
│   ├── §3  print goal + box info
│   ├── §4  plot "Original" tab (mechanism @ θ2=0 + Pin A path)
│   ├── §5  plot "Original Angles" tab (8 frames every 15°)
│   ├── §6  write per-5° box table to Excel (Original sheet)
│   ├── §7  design search: triple nested loop over (theta1, dP1A, dP2A)
│   ├── §8  sort results by area, plot Best, write Best sheet
│   └── §9  fallback path: if no candidate satisfies constraints,
│           plot two hand-picked alternates
│
└── local functions
    ├── solveDesign(...)        — Newton-Raphson on 6 loop equations, all θ2
    ├── boxFromPts(pts)         — axis-aligned bounding box for a set of points
    ├── plotMechanism(...)      — draws skeleton + bounding box for one frame
    ├── plotAnglesFigure(...)   — tiled multi-frame mechanism plot
    └── printAngleTable(...)    — per-5° table to console + Excel sheet
```

A single tabbed `uitabgroup` figure holds all four result tabs:
`Original`, `Original Angles`, `Best`, `Best Angles`. After all plots are drawn,
every mechanism axes is forced to share `xlim`/`ylim` (computed from the union of
joint positions across both designs and all `θ2`) so the original and best designs
are visually comparable, then the figures are re-exported.

---

## 4. The kinematics solver (`solveDesign`)

This is the core. It loops `theta2` from `0°` to `120°` in `1°` steps and runs
Newton-Raphson at each step to find the six unknown angles:

```
x = [theta23, theta14, theta46, theta36, theta8, theta7]
```

The residual `f(x)` is built from **three vector loops** of the linkage (six scalar
equations — x and y components of each loop). The Jacobian `J` is the 6×6 matrix
of partials of `f` w.r.t. `x`. Each NR step is the standard

```
dx     = J \ f
x_new  = x - dx
x_new  = atan2(sin x_new, cos x_new)   % wrap into (-pi, pi]
```

with convergence test `||f||_inf < 1e-8` AND `||dx||_inf < 1e-8`, capped at 200
iterations. The angle wrap keeps NR from drifting onto a different branch.

**Initial guess strategy.** This is the subtle bit:

- For the original design, `x0` is hard-coded to a guess close to the known
  configuration at `θ2 = 0°` (line 56). Without it, NR can converge to a different
  posture.
- For every subsequent `θ2`, `x` carries over from the previous step — the
  previous solution is the initial guess for the next, which keeps NR on the same
  branch as `θ2` advances.
- For the design-space search, the initial guess for **every** candidate is
  `x0_base = allAngles.raw(:,1)` (the converged angles at `θ2 = 0°` from the
  original design). Reusing the same warm start across candidates is what makes
  the search practical.

After NR converges, the joint positions are reconstructed by walking the linkage
from `O2` outward (lines 400–409). Pin A is then placed by:

```
u = (P2 - P1) / R8                 % unit vector along P1->P2
v = [-u(2), u(1)]                  % CCW perpendicular
PinA = P1 + pinA_xL*u + pinA_yL*v  % local-frame coords from law of cosines
```

where `pinA_xL` and `pinA_yL` come from the rigid triangle (P1, P2, PinA) sides
`R8`, `dP1A`, `dP2A`. `pinA_yL` carries the `pinASide` sign.

Returns:
- `allAngles.raw` — 6×N matrix of converged angles (N = 121)
- `pinA.Ax`, `pinA.Ay` — Pin A trajectory
- `boxInfo` — bounding box AT `θ2 = 0`, plus `joints` (snapshot at θ2=0),
  `labels`, and `allJoints` (11×2×121 tensor of every joint at every θ2)

---

## 5. The design search (the big nested loop)

```
for iT  = theta1_range          % 21 values, ±15° around theta1_orig
  for iD1 = dP1A_range          % 16 values, 0.8..1.1 × dP1A_orig
    for iD2 = dP2A_range        % 17 values, 0.7..1.1 × dP2A_orig
        skip if triangle inequality fails for (R8, dP1A, dP2A)
        try
            solveDesign(...)
            if minDist(PinA, goal) < goalTol  AND  area < origArea
                store the candidate
            end
        catch
            % NR did not converge — silently skip
        end
    end
  end
end
```

That's `21 × 16 × 17 = 5712` candidates. For each one we re-solve the full
`θ2 = 0..120°` sweep. The two filters that decide whether a candidate is kept:

1. **Goal reachability:** `min |PinA(θ2) − goal|` over the sweep must be below
   `goalTol = 15 mm`.
2. **Box improvement:** the bounding box area at `θ2 = 0°` must be strictly
   smaller than the original.

Surviving candidates are sorted ascending by `area` and `results(1)` is the best.

---

## 6. Outputs

Files dropped next to the script (not in `notes/`):

- `ENME473_Q7_Original.png` — original mechanism @ θ2=0 + Pin A path
- `ENME473_Q7_Original_Angles.png` — original at θ2 = 0:15:120
- `ENME473_Q7_BestDesign.png` — best mechanism @ θ2=0 + Pin A path
- `ENME473_Q7_BestDesign_Angles.png` — best at θ2 = 0:15:120
- `ENME473_Q7_BoxDimensions.xlsx` — two sheets (`Original`, `Best`), each with
  `Theta2_deg, Width_mm, Height_mm, Area_mm2` every 5° from 0 to 120

Console output: original design summary, search progress, best/alternate design
summary, and a comparison table.

---

## 7. Reproducing this from scratch

If you wanted to write the same script for a different mechanism, the recipe is:

1. **Inventory** every link length and offset angle of the mechanism. Decide which
   are fixed and which become design parameters.
2. **Write the loop equations** by hand. For an N-loop linkage with M unknown
   angles you need N vector equations → 2N scalar equations, and you should have
   M = 2N.
3. **Symbolically (or by hand) differentiate** the residual vector to get the
   Jacobian. Drop both into a per-step Newton-Raphson with angle wrapping.
4. **Find one good initial guess** at the starting input angle (sketch it; measure
   off the figure; whatever). Then carry the converged solution forward as the
   guess for the next step — that's how you stay on a single branch as the input
   sweeps.
5. **Reconstruct joint positions** by walking the linkage from a fixed pivot
   outward, using the now-known angles.
6. **Define what "good" means** (here: small bounding box, reachable goal). Wrap
   the solver in a function that returns a metric.
7. **Sweep the design parameters** with nested `for` loops, guarding each call in
   `try/catch` to silently skip non-convergent candidates and any that violate
   geometric feasibility (e.g. the triangle inequality on the link-8 triangle).
8. **Sort, plot, export.** Tabbed figure for interactive review; PNGs and an
   `xlsx` workbook for the report.

---

## 8. Things to know if you tweak this

- **`pinASide = ±1`** flips Pin A across the P1–P2 axis. If your plotted mechanism
  is mirrored relative to the assignment figure, flip this first.
- **`x0` (line 56)** is fragile — change the geometry too much and NR may converge
  to a different posture or fail. If the search returns zero candidates, sanity
  check that the original solve still produces a sensible plot.
- **`goalTol = 15 mm`** is the only fitness-relaxation knob. Tighten it for a
  stricter match; loosen if the search returns nothing.
- The search ranges (`linspace(...)`) are deliberately asymmetric on `dP1A`/`dP2A`
  because shrinking these tends to reduce the box more than expanding them — feel
  free to widen them if you have CPU to burn.
- The nested loop is embarrassingly parallel; swap the outer `for` for `parfor`
  if you want to cut runtime, but you'll need to handle the `results = [results; r]`
  growth pattern (preallocate or use a sliced cell).

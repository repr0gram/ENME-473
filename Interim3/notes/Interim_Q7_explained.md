# Interim_Q7.m — Walkthrough

This note explains `Interim3/Interim_Q7.m` end-to-end so someone unfamiliar with
the script can reproduce the workflow: solve the original 8-link mechanism, run
an **inverse** design search that places Pin A exactly on the goal by
construction, and produce plots + an Excel report.

---

## 1. What the script is doing (high level)

The mechanism is an 8-link planar linkage driven by an input crank `theta2`. We
care about the trajectory of **Pin A** — a point fixed on the rigid triangle of
link 8. The goal is:

1. Solve the kinematics for `theta2 = 0..120°` with the original geometry.
2. Record the **goal**: where Pin A sits at `theta2 = 120°`.
3. Compute the **enclosing (axis-aligned) bounding box** of the whole mechanism
   at `theta2 = 0°`.
4. Find an alternative configuration whose mechanism has a *smaller* box at
   `theta2 = 0°` while Pin A still passes through the goal exactly.

The three design knobs are:
- `theta1` — angle of the ground link (`|R1|` is held fixed; `O2` stays at the origin)
- `dP1A` — distance from Pin A to the *far* link-8 pin (P1)
- `dP2A` — distance from Pin A to the *near* link-8 pin (P2)

`R8` (the distance between the two link-8 pins themselves) is fixed.

### The key idea: work backwards from the goal

A naïve forward search sweeps all three design parameters on a grid and
asks "did Pin A happen to pass near the goal?" for each candidate. That
wastes most of the compute on designs that fail the goal check.

The inverse approach exploits a structural property of the linkage: in the
Newton-Raphson loop-equation solver, **`dP1A` and `dP2A` do not appear in the
posture equations at all**. They only show up afterward, when Pin A is placed on
the already-solved P1 and P2 pins of link 8. That means, for a fixed `theta1`,
the joint positions of every link except the Pin A placement are **independent**
of the link-8 triangle shape.

So the script does this instead:

1. For each candidate `theta1` and each target input angle `theta2_target`, solve
   the mechanism. Extract `P1` and `P2` at that posture.
2. The goal is a fixed point; the required triangle side lengths follow directly
   from two distance queries:
   ```
   dP1A_required = |P1 - goal|
   dP2A_required = |P2 - goal|
   ```
3. Determine which side of the `P1→P2` axis the goal lies on (that's
   `pinASide`).
4. Re-place Pin A on the joint tensor with these triangle sides and evaluate
   the bounding box at `theta2 = 0°`.

Every candidate that survives the feasibility checks lands Pin A *exactly* on
the goal (verified numerically at ~7 × 10⁻¹² mm residual). Only `theta1` and
`theta2_target` are swept — `dP1A` and `dP2A` are *outputs* of the goal
constraint, not inputs to search.

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

Pin A is **off** the P1–P2 line; its position on link 8 is determined by the
rigid triangle with sides `R8`, `dP1A`, `dP2A`. The sign `pinASide = ±1` picks
which side of the P1→P2 axis Pin A sits on.

---

## 3. File / code structure

```
Interim_Q7.m
├── (script body)
│   ├── §1  fixed link lengths + nominal angles (alpha, beta, gamma, theta1)
│   ├── §2  solve original design (1 call to solveDesign)
│   ├── §3  print goal + box info, set origArea baseline
│   ├── §4  plot "Original" tab (mechanism @ θ2=0 + Pin A path)
│   ├── §5  plot "Original Angles" tab (8 frames every 15°)
│   ├── §6  write per-5° box table to Excel (Original sheet)
│   ├── §7  INVERSE SEARCH (2D grid: theta1 × theta2_target)
│   ├── §8  sort by area, re-solve + plot Best, write Best sheet
│   └── §9  report alternate design
│
└── local functions
    ├── replacePinA(...)       — re-place Pin A on an existing joint tensor
    │                            using a new (dP1A, dP2A, pinASide) triangle
    ├── solveDesign(...)       — Newton-Raphson on 6 loop equations, all θ2
    ├── boxFromPts(pts)        — axis-aligned bounding box for a set of points
    ├── plotMechanism(...)     — draws skeleton + bounding box for one frame
    ├── plotAnglesFigure(...)  — tiled multi-frame mechanism plot
    └── printAngleTable(...)   — per-5° table to console + Excel sheet
```

A single tabbed `uitabgroup` figure holds all four result tabs:
`Original`, `Original Angles`, `Best`, `Best Angles`. After all plots are drawn,
every mechanism axes is forced to share `xlim`/`ylim` (computed from the union of
joint positions across both designs and all `θ2`) so the original and best
designs are visually comparable, then the figures are re-exported.

---

## 4. The kinematics solver (`solveDesign`)

This is the core. It loops `theta2` from `0°` to `120°` in `1°` steps and runs
Newton-Raphson at each step to find the six unknown angles:

```
x = [theta23, theta14, theta46, theta36, theta8, theta7]
```

The residual `f(x)` is built from **three vector loops** of the linkage (six
scalar equations — x and y components of each loop). The Jacobian `J` is the 6×6
matrix of partials of `f` w.r.t. `x`. Each NR step is the standard

```
dx     = J \ f
x_new  = x - dx
x_new  = atan2(sin x_new, cos x_new)   % wrap into (-pi, pi]
```

with convergence test `||f||_inf < 1e-8` AND `||dx||_inf < 1e-8`, capped at 200
iterations. The angle wrap keeps NR from drifting onto a different branch.

**Important observation.** Look at the six residual equations inside
`solveDesign`. None of them involve `dP1A`, `dP2A`, or `pinASide`. They depend
only on the link lengths `R1..R8, R14, R23, R36, R46`, the offset angles
`alpha, beta, gamma`, and the input `theta2`. This is what makes the inverse
search possible: once `theta1` is fixed, you can solve for every non-PinA joint
*without* knowing the link-8 triangle shape, then figure out the triangle
afterward.

**Initial guess strategy.** This is the subtle bit:

- For the original design, `x0` is hard-coded to a guess close to the known
  configuration at `θ2 = 0°`. Without it, NR can converge to a different
  posture.
- For every subsequent `θ2`, `x` carries over from the previous step — the
  previous solution is the initial guess for the next, which keeps NR on the
  same branch as `θ2` advances.
- For the inverse search, the initial guess for **every** `theta1` candidate is
  `x0_base = allAngles.raw(:, 1)` (the converged angles at `θ2 = 0°` from the
  original design). Reusing the same warm start across candidates is what makes
  the search fast.

After NR converges, the joint positions are reconstructed by walking the
linkage from `O2` outward. Pin A is then placed by:

```
u = (P2 - P1) / R8                 % unit vector along P1->P2
v = [-u(2), u(1)]                  % CCW perpendicular
PinA = P1 + pinA_xL*u + pinA_yL*v  % local-frame coords from law of cosines
```

where `pinA_xL` and `pinA_yL` come from the rigid triangle with sides `R8`,
`dP1A`, `dP2A`. `pinA_yL` carries the `pinASide` sign.

Returns:
- `allAngles.raw` — 6×N matrix of converged angles (N = 121)
- `pinA.Ax`, `pinA.Ay` — Pin A trajectory
- `boxInfo` — bounding box AT `θ2 = 0`, plus `joints` (snapshot at θ2=0),
  `labels`, and `allJoints` (11×2×121 tensor of every joint at every θ2)

---

## 5. The inverse design search

```
for iT = theta1_range                  % 81 values, ±20° around theta1_orig
    solveDesign(...) once with any dummy (dP1A, dP2A)
    jointsRaw = allJoints              % mechanism shape, independent of dP1A/dP2A

    for iA = theta2_target_range       % 31 values, 90..120 deg
        idx = theta2_target + 1
        P1 = jointsRaw(8, :, idx)
        P2 = jointsRaw(9, :, idx)

        dP1A_req = |goal - P1|
        dP2A_req = |goal - P2|

        skip if dP1A_req or dP2A_req outside [0.5, 1.5] × original
        skip if triangle inequality fails for (R8, dP1A_req, dP2A_req)

        compute pinASide from sign(dot(goal - P1, perp(P2 - P1)))
        replacePinA(jointsRaw, dP1A_req, dP2A_req, pinASide, R8)
        compute bounding box at theta2 = 0 from the updated joint tensor

        if area < origArea
            store the candidate
        end
    end
end
```

That's `81 × 31 = 2511` candidates. For each `theta1` we run the Newton-Raphson
kinematics solver once; the inner loop over `theta2_target` does only
arithmetic on already-solved joint positions. The filters that decide whether a
candidate is kept:

1. **Physical bounds on the triangle sides:** `0.5 × dP1A_orig ≤ dP1A_req ≤ 1.5 × dP1A_orig`
   (same for `dP2A`). This rejects ridiculously large or small link-8 triangles.
2. **Triangle inequality:** `dP1A_req + dP2A_req > R8` and
   `|dP1A_req - dP2A_req| < R8`, so the rigid triangle actually closes.
3. **Box improvement:** bounding box area at `θ2 = 0°` must be strictly smaller
   than the original.

Surviving candidates are sorted ascending by area and `results(1)` is the best.

Note that **there is no goal-reach tolerance** — the required triangle sides
are computed directly from the goal location, so Pin A lands on the goal by
construction. The script verifies this at the end by re-solving the best
design through the normal `solveDesign` path and checking the residual
(~7 × 10⁻¹² mm in practice, i.e. machine precision).

---

## 6. The `replacePinA` helper

After the joint tensor is computed for a given `theta1`, changing the link-8
triangle shape only affects Pin A (row 11). `replacePinA` walks every `theta2`
frame, reads `P1` and `P2` from rows 8 and 9, and computes Pin A from the new
`(dP1A, dP2A, pinASide)`:

```matlab
pinA_xL = (dP1A^2 - dP2A^2 + R8^2) / (2*R8);    % law of cosines
pinA_yL = pinASide * sqrt(max(0, dP1A^2 - pinA_xL^2));
for k = 1:N
    u = (P2(k) - P1(k)) / R8;
    v = [-u(2), u(1)];
    joints(11, :, k) = P1(k) + pinA_xL*u + pinA_yL*v;
end
```

No NR, no new solve — just a fixed transform applied to the joint array.

---

## 7. Outputs

Files dropped next to the script:

- `ENME473_Q7_Original.png` — original mechanism @ θ2=0 + Pin A path
- `ENME473_Q7_Original_Angles.png` — original at θ2 = 0:15:120
- `ENME473_Q7_BestDesign.png` — best mechanism @ θ2=0 + Pin A path
- `ENME473_Q7_BestDesign_Angles.png` — best at θ2 = 0:15:120
- `ENME473_Q7_BoxDimensions.xlsx` — two sheets (`Original`, `Best`), each with
  `Theta2_deg, Width_mm, Height_mm, Area_mm2` every 5° from 0 to 120

Console output: original design summary, search progress, best/alternate design
summary with residual error, and a comparison table.

---

## 8. Reproducing this from scratch

If you wanted to write the same script for a different mechanism, the recipe is:

1. **Inventory** every link length and offset angle. Decide which are fixed and
   which become design parameters.
2. **Write the loop equations** by hand. For an N-loop linkage with M unknown
   angles you need N vector equations → 2N scalar equations, and you should
   have M = 2N.
3. **Symbolically (or by hand) differentiate** the residual vector to get the
   Jacobian. Drop both into a per-step Newton-Raphson with angle wrapping.
4. **Identify which design parameters show up in the loop equations and which
   don't.** Any parameter that only affects a terminal point (like Pin A on a
   rigid triangle) can be inverted after the fact, exactly like the Pin A
   triangle here. Any parameter that *does* appear in the loop equations has to
   stay in the outer search.
5. **Find one good initial guess** at the starting input angle. Then carry the
   converged solution forward as the guess for the next step — that's how you
   stay on a single branch as the input sweeps.
6. **Reconstruct joint positions** by walking the linkage from a fixed pivot
   outward, using the now-known angles.
7. **Define the objective** (here: minimize bounding box at `θ2 = 0°`). Put it
   inside a helper that runs on a joint tensor.
8. **Sweep only the parameters that touch the loop equations.** For every
   candidate, solve the mechanism once, then do a cheap inner loop over any
   "attach point" parameters that you can compute directly from the goal.
9. **Feasibility prefilters.** Triangle inequality for rigid triangles, sign
   conventions for attach-side choices, and `try/catch` around the NR solver to
   silently drop non-convergent candidates.
10. **Sort, plot, export.** Tabbed figure for interactive review; PNGs and an
    `xlsx` workbook for the report.

---

## 9. Things to know if you tweak this

- **`x0` (the hard-coded NR start)** is fragile — change the geometry too much
  and NR may converge to a different posture or fail. If the search returns
  zero candidates, sanity check that the original solve still produces a
  sensible plot.
- **`theta1_range` and `theta2_target_range`** are the only two sweep knobs.
  Widen them for more exploration; 81 × 31 = 2511 candidates takes a few
  seconds on a laptop and is plenty dense in practice.
- **`dP1A_bounds`, `dP2A_bounds`** reject designs that would require a very
  different link-8 triangle than the original. Loosen them if you want to
  explore more aggressive geometries.
- **There is no goal tolerance** — if you want a softer constraint (e.g., "the
  mechanism is allowed to stop slightly short of the goal"), you'd have to
  reintroduce a `goalTol` check on the forward distance. The inverse approach
  has no equivalent knob because it solves exactly.
- The outer `for iT` loop is embarrassingly parallel; swap it for `parfor` if
  you want to cut runtime, but you'll need to handle the
  `results = [results; r]` growth pattern (preallocate or use a sliced cell).

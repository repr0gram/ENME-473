# 7 Design Optimization

Question 7 asks us to treat the mechanism as a starting point for a simple
design study. With the input crank at 120°, Pin A reaches a specific "goal"
location; with the input crank at 0°, the whole mechanism is bounded by a
rectangular envelope whose area we would like to make as small as possible. The
task is to vary three design parameters — the location of one ground pivot and
the two link-8 triangle side lengths that fix Pin A relative to the rest of
link 8 — and find an alternative geometry that reaches the goal while occupying
a smaller enclosing box at the home position.

Rather than sweeping the three design parameters on a blind grid and checking
after the fact whether Pin A happened to pass near the goal, the approach taken
here inverts the problem: the goal position is treated as a constraint, and the
two link-8 triangle side lengths are solved for analytically so that Pin A
lands on the goal by construction. Only the ground-link angle and the target
input angle remain as free variables. Every candidate in the search therefore
reaches the goal *exactly* — no tolerance is needed — and the entire search
budget is spent on minimizing the enclosing box.

The investigation below sets up the problem, explains the inverse search,
reports the original and best-of-search geometries, compares them side by side,
and discusses the practical implications of the result.

## 7.1 Problem Setup and Design Parameters

The nominal mechanism is the same 8-link planar linkage analyzed in Parts I
and II. All link lengths and ternary offset angles (α = 3.2°, β = 3.5°,
γ = 1.8°) are held fixed, as is the length of the ground link (|R₁|) and the
length R₈ between the two pins on link 8. Pin A is a rigid point on link 8 that
does **not** lie on the line joining those two pins; its position is fully
determined once the two triangle sides dP1A and dP2A are chosen.

The **goal** is defined as the location of Pin A when the input angle
θ₂ = 120°, computed from the original design:

- Goal coordinates: (x, y) = (−1580.75, 552.07) mm

The **enclosing box** is the axis-aligned rectangle that contains every joint
of the mechanism at θ₂ = 0°. Its width, height and area for the original
design are reported in §7.3.

Three design parameters are allowed to vary:

| Symbol | Meaning                                              | Original value |
|--------|------------------------------------------------------|----------------|
| θ₁     | Ground-link angle (\|R₁\| fixed, O₂ stays at origin) | 163.42°        |
| dP1A   | Distance from Pin A to the far link-8 pin (P1)       | 580.70 mm      |
| dP2A   | Distance from Pin A to the near link-8 pin (P2)      | 336.90 mm      |

Varying θ₁ rotates the ground pivot O₄ about O₂ along a circle of radius |R₁|,
so O₄ moves but the ground-link length does not change. Varying dP1A and dP2A
reshapes the rigid P1–P2–PinA triangle on link 8; R₈ is kept constant, so this
amounts to sliding Pin A around on link 8 without altering where the two
link-8 pins sit relative to each other.

## 7.2 Methodology

The design search exploits a structural observation about the linkage: the
Newton–Raphson loop-equation solver used for the posture analysis depends only
on link lengths and the ground-link angle θ₁ — it does **not** depend on
dP1A or dP2A. Those two parameters are applied only *afterward*, when Pin A is
placed on the already-solved link-8 pins P1 and P2 via the rigid P1–P2–PinA
triangle. For a fixed θ₁, the positions of P1 and P2 at every input angle are
therefore independent of how link 8 is shaped.

This lets us flip the problem around. Rather than guessing (dP1A, dP2A) and
checking whether Pin A ends up near the goal, we **require** Pin A to land on
the goal at some input angle θ₂_target and solve for the triangle sides
directly:

    dP1A_required  =  |P1(θ₂_target) − goal|
    dP2A_required  =  |P2(θ₂_target) − goal|

Whichever side of the P1→P2 axis the goal lies on determines the sign of the
third triangle coordinate (pinASide = ±1). The resulting triangle exactly
reproduces Pin A = goal at θ₂ = θ₂_target, to machine precision.

**Search grid.** Only two parameters are swept:

- θ₁: 81 values in [θ₁₀ − 20°, θ₁₀ + 20°] (every 0.5°)
- θ₂_target: 31 integer values in [90°, 120°]

The total is 81 × 31 = 2511 candidate designs. For each candidate the full
mechanism is solved once (the Newton–Raphson sweep over θ₂ = 0°–120° in 1°
steps) to obtain the joint history, then the inner loop over θ₂_target reduces
to distance queries on already-computed positions.

**Feasibility filters.** Each candidate must satisfy three conditions:

1. The computed dP1A and dP2A must lie within ±50% of the original values (a
   loose physical-plausibility bound on the link-8 triangle).
2. The triangle inequality must hold on (R₈, dP1A, dP2A) so that the rigid
   triangle closes.
3. The bounding box area at θ₂ = 0° must be strictly smaller than the
   original. The box is computed from the full joint set after Pin A has been
   re-placed with the new triangle sides.

**Reporting.** All surviving candidates are sorted ascending by box area at
θ₂ = 0°; the top result is reported as the "best" design and the runner-up as
an alternate. For both the original and the best design, the enclosing box
width, height and area are tabulated every 5° from θ₂ = 0° to 120° and
exported to `ENME473_Q7_BoxDimensions.xlsx`. Mechanism snapshots at 15°
increments are plotted side-by-side on shared axes so that the envelope of
motion can be compared visually.

A sanity check is performed on the best design: the mechanism is re-solved
through the normal solver with the computed (dP1A, dP2A, pinASide), and the
distance between Pin A at θ₂ = θ₂_target and the original goal is measured.
For the reported best design this residual is **7.08 × 10⁻¹² mm**, i.e. the
goal is hit to machine precision.

## 7.3 Results

**Original design.** At θ₂ = 0° the mechanism occupies a 721.18 × 317.25 mm
rectangle with an area of 228,794.84 mm². As the input crank rotates, the
envelope grows — the box reaches its maximum area of roughly 1.12 × 10⁶ mm²
near θ₂ = 100°, driven almost entirely by the horizontal extent of Pin A's
trajectory sweeping to the left (see `ENME473_Q7_Original_Angles.png`). At
θ₂ = 120° Pin A is at (−1580.75, 552.07) mm; this position is fixed as the
optimization goal.

**Design search.** Of the 2511 candidate (θ₁, θ₂_target) combinations
evaluated, **177** produced a valid design that (i) satisfies the physical
triangle bounds, (ii) closes the rigid triangle, and (iii) has a smaller
enclosing box at θ₂ = 0° than the original. All 177 reach the goal exactly.

**Best design.** Sorting by box area, the best design is:

| Parameter | Value               | Original           | Change      |
|-----------|---------------------|--------------------|-------------|
| θ₁        | 165.42°             | 163.42°            | +2.00°      |
| O₄ (x, y) | (−239.52, 62.28) mm | (−237.20, 70.60) mm | −8.32 mm in y |
| dP1A      | 668.51 mm           | 580.70 mm          | +87.81 mm   |
| dP2A      | 414.20 mm           | 336.90 mm          | +77.30 mm   |

The enclosing box at θ₂ = 0° shrinks from 721.18 × 317.25 mm (228,795 mm²) to
**832.47 × 137.12 mm (114,147 mm²)**, a **50.11%** reduction in area. The width
of the new box is larger (+111.29 mm), but the height drops by more than half
(from 317.25 mm to 137.12 mm), and the height reduction dominates. Pin A on
the best design passes through the goal exactly at θ₂ = **112°**, with a
residual error of 7.08 × 10⁻¹² mm.

**Alternate design.** The second-best candidate in the sorted list keeps the
same ground pivot as the best design and differs only in the link-8 triangle
and the target input angle:

| Parameter | Value               |
|-----------|---------------------|
| θ₁        | 165.42°             |
| O₄ (x, y) | (−239.52, 62.28) mm |
| dP1A      | 680.70 mm           |
| dP2A      | 426.30 mm           |
| Box       | 844.50 × 137.12 mm  |
| Area      | 115,797 mm²         |
| Reduction | 49.39%              |
| Goal pass | θ₂ = 111° (exact)   |

The alternate demonstrates that the optimum sits inside a family of closely
related designs, not at an isolated numerical artefact.

**Figures.** `ENME473_Q7_Original.png` and `ENME473_Q7_BestDesign.png` show
the θ₂ = 0° posture of each design with its enclosing box and Pin A's full
trajectory superimposed. `ENME473_Q7_Original_Angles.png` and
`ENME473_Q7_BestDesign_Angles.png` show both mechanisms at 15° increments
across the full input range, on shared axes, so that the evolution of the
envelope can be compared frame by frame.

## 7.4 Comparison

The table below summarizes the three designs side by side.

| Design    | θ₁ (°)  | O₄ x (mm) | O₄ y (mm) | dP1A (mm) | dP2A (mm) | Box W × H (mm)     | Box Area (mm²) | Area Δ   | Goal @ θ₂ |
|-----------|---------|-----------|-----------|-----------|-----------|--------------------|----------------|----------|-----------|
| Original  | 163.42  | −237.20   | 70.60     | 580.70    | 336.90    | 721.18 × 317.25    | 228,795        | —        | 120°      |
| Best      | 165.42  | −239.52   | 62.28     | 668.51    | 414.20    | 832.47 × 137.12    | 114,147        | −50.11% | 112°      |
| Alternate | 165.42  | −239.52   | 62.28     | 680.70    | 426.30    | 844.50 × 137.12    | 115,797        | −49.39% | 111°      |

A few observations:

- Both optimized designs grow wider than the original (by roughly 110 mm), but
  shrink in height by a much larger amount (from 317 mm to 137 mm). The
  dominant vertical extent of the original posture at θ₂ = 0° comes from
  link 8 carrying Pin A far above the P1–P2 line. Enlarging dP1A and dP2A
  *together* rotates the Pin A attachment on link 8 toward a flatter
  configuration at the home position, bringing the top of the mechanism much
  closer to the P1–P2 line. Roughly 180 mm of height reduction for 110 mm of
  width increase is a very favorable trade against an axis-aligned box.
- The ground-pivot shift (+2° in θ₁, about 8 mm downward at O₄) is modest.
  Most of the improvement comes from reshaping the link-8 triangle, not from
  moving the frame.
- Neither optimized design reaches the goal at exactly 120°; both reach it
  slightly earlier (at 112° and 111° respectively). This is allowed by the
  project statement, which explicitly notes that the mechanism need not reach
  the goal at an input angle of exactly 120°.
- Across the full θ₂ range, the best design's peak box area
  (≈ 1.052 × 10⁶ mm² near θ₂ = 100°) is slightly smaller than the original's
  peak (≈ 1.123 × 10⁶ mm² near θ₂ = 100°), so the improvement is not confined
  to the home position — the optimized mechanism is modestly more compact
  throughout the motion. The largest gain, however, is very clearly at
  θ₂ = 0°, which is the objective.

## 7.5 Discussion

The optimization meets both goals set out in the project statement: reach the
original goal position exactly and reduce the enclosing-box area at θ₂ = 0°.
A 50% reduction in the home-position envelope area is a substantial result and
qualifies for the bonus described in the Q7 prompt.

A few points about the method and its limits are worth noting.

**Why the inverse formulation is exact.** Pin A's coordinates depend on
(dP1A, dP2A, pinASide) only through the rigid link-8 triangle; the rest of the
mechanism is unaffected. This means that once θ₁ and the target input angle
are fixed, the required triangle sides are determined by two distance queries
against the goal — no iteration, no tolerance, no approximation. The only
unknowns left in the search are θ₁ and the target angle, so the 2D grid is
both much cheaper and much more successful than a 3D grid over
(θ₁, dP1A, dP2A) would be.

**Search density and bounds.** The grid is 81 × 31 = 2511 candidates. A
denser grid or a wider range would find marginally better designs, but 177
valid candidates already cluster tightly around the reported best, which
suggests the optimum is well captured. The triangle-side bounds of ±50% of the
original values were chosen to reject geometries that would require
substantially re-building link 8; loosening them would likely yield further
improvement at the cost of more dramatic hardware changes.

**Objective scope.** Only the enclosing box at θ₂ = 0° is used as the
objective. The full-range data in `ENME473_Q7_BoxDimensions.xlsx` shows that
the envelope area grows with θ₂, peaking near 100° for both designs; if the
relevant constraint for a given application were the worst-case envelope
across the full cycle rather than at the home position, the optimal design
could look different. For this problem, the home-position box is what the
Q7 prompt asks for, and that is what is minimized here.

**What is not checked.** Kinematic feasibility (closure and triangle
inequality) is the only constraint enforced. No consideration is given to
interference between links, transmission-angle degradation, or the effect of
the modified geometry on the mechanism's velocity and acceleration profiles.
A next step would be to re-run Parts I and II on the best design to verify
that its dynamic behavior remains acceptable. A stricter check on how far the
link-8 triangle has been reshaped (dP1A is about 15% longer and dP2A about 23%
longer than their nominal values) would also be in order before committing to
the new geometry.

Overall, the inverse formulation turned a moderately successful approximate
search into an exact one while simultaneously cutting the computational cost,
and the resulting best design halves the enclosing-box area at the home
position — a significantly larger improvement than a forward grid search
would have found in the same parameter ranges.

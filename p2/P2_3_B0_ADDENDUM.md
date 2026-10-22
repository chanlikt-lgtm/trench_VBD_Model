# Addendum — B0 gate and B0→B1 area (closes the pending items)

**Date:** 2026-10-22
**Applies to:** `P2_3_SESSION_B0_OSQP_BODY_CUT.md` (filed verbatim, unedited).
**Full results record:** `P2_3_B0_B1_BAND.md`; commit `438d7cc`.

The session record was written while B0 was still polishing. Two of its
pending items are now closed.

## 1. B0 full-set F7 — PASS

```
active rows        4,435
stationarity       5.824e-11
complementarity    5.016e-17
min slack         -8.147e-15
violations         0
lambda min        +0.000e+00
oxide residual     4.083e-13
overshoot         +4.268e-15
GATE              PASS
```

No active-set repair was required (unlike B1, where two negative multipliers
were dropped).

## 2. Bounded neutral area and the B0→B1 increment

```
B0  [0.021644935001, 0.021649266844]
B1  [0.021987652318, 0.021991846424]

B0 -> B1   [+3.383855e-04, +3.469114e-04]   +1.563 % … +1.603 %
```

Strictly positive, non-overlapping, and 1.5× the entire L1b→L1 step.

## 3. Applying the record's own decision gate

The record specified: *large, non-overlapping positive increment → strong
response to the band intervention; prioritise exact donor integration before
assigning the magnitude solely to boundary resolution.*

Measured against that gate:

- area change: 1.56–1.60 % relative;
- donor change: 1.54e-4 relative (ΔΣf = +0.3035).

The two differ by a factor of ~100, so the centroid-source term cannot account
for the magnitude. It remains a real confound at the 0.015 % level, and the
exact partial-element load remains the control that would remove it — but it is
not required to conclude that boundary resolution drives this increment.

The residual question the gate leaves open is narrower: what *fraction* of the
1.58 % is source-driven. Exact integration answers it; it does not change the
sign or the order of magnitude.

## 4. Front cuts (ε = 1e-8)

```
x = 2.35   B0  (-0.0649, +0.0685)      B1  (-0.0746,-0.0727)+(-0.0698,+0.0757)
x = 2.41   B0  (-0.2026, +0.1957)      B1  (-0.2017,-0.2005)+(-0.1972,+0.1975)
```

The band widens the neutral region upstream by 0.0097 µm at x = 2.35, consistent
with the area increase. B1's satellite spans remain unverified as physical
features; they sit at the band edge and merge by ε = 1e-5.

## 5. Status

B0 certified; B1 certified; B2 not authorised; P2.3 open; P1 frozen; Stage 3
suspended. No convergence order is claimed: B0→B1 is a single jump of factor 2.1
in h.

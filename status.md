# Lead--Lag Control Lab -- Development Status

## Repository State (Feb 2026)

-   Stable version is currently deployed via `gh-pages` branch.
-   Development (challenge version) is on `main` branch locally.
-   Challenge system implemented but not fully completed.

------------------------------------------------------------------------

# Implemented

## Core Lab (Stable)

-   Plant G(s) input (num/den)
-   Lead--lag controller
-   Bode (before/after)
-   Nyquist (before/after)
-   Time simulation (RK4)
-   Adjustable timestep
-   Axis scaling (frequency ×10\^k, time scaled inversely)
-   Closed-loop pole computation (Durand--Kerner)
-   Plotly + numeric.js bundled locally
-   GitHub Pages deployment working

## Challenge System (Phase 1.5--2)

### Infrastructure

-   Challenge selector
-   Margin targets (Am, Phim)
-   Score logic (0 unless targets met)
-   Status panel with per-metric feedback
-   Time-domain mask overlay on output plot
-   Mask violation detection
-   Mask integrated into score gating
-   Infinite gain margin case handled correctly

### Implemented Challenges

-   ✅ Challenge 1: Robust & Calm (should change name)
    -   Margin targets active
    -   Time-domain mask implemented:
        -   Overshoot limit (y \< 1.2 for 0 \< t \< 5)
        -   Settling band 0.95--1.05 for 5 \< t \< 20
        -   Load disturbance band 0.7--1.3 for 20 \< t \< 40
-   ⚠ Challenge 2: Fast but Safe (change)
    -   Structure exists
    -   Not yet fully defined (mask + tuned targets missing)
-   ⚠ Challenge 3: Disturbance Fighter (change)
    -   Structure exists
    -   Not yet fully defined (disturbance metrics + mask missing)

------------------------------------------------------------------------

# Known Limitations / TODO

## 1. Responsive Display Scaling

-   Plot sizes currently fixed.
-   Should adapt to:
    -   Screen resolution
    -   Laptop vs large monitor
    -   Possibly mobile (optional)
-   Consider dynamic height based on window.innerHeight.

## 2. Measurement Noise Metric

Idea: - Compute RMS of measurement noise contribution at output. - Add
as optional performance penalty.

Considerations: - Must avoid heavy computation. - Could reuse simulated
time data. - RMS(y_noise_component) over specified window.

Alternative (if too expensive): - Penalize high controller gain at high
frequency. - Evaluate \|C(jω)\| above some ω threshold. - Add soft
penalty if gain too large.

Educational tradeoff: - Keep conceptual clarity. - Avoid black-box
scoring.

## 3. Challenge 2 Design Ideas

-   Tight settling time constraint.
-   Limit overshoot strictly.
-   Possibly penalize excessive control effort.

## 4. Challenge 3 Design Ideas

-   Disturbance rejection metric.
-   Evaluate peak deviation after load step.
-   Evaluate recovery time to ±5% band.

## 5. UI Refinements

-   Improve visual hierarchy of status panel.
-   Consider color coding (green/red).
-   Consider highlighting active constraints more clearly.
-   Make **score** stand out more, visually
-   Remove some text clutter 

------------------------------------------------------------------------

# Architectural Notes

-   Manual state handling (no framework).
-   Plotly.react used consistently.
-   Frequency/time scaling tightly coupled.
-   Mask system robust and generalizable.
-   Gain margin infinity case now handled properly.

------------------------------------------------------------------------

# Future Publishing Workflow

Stable deployment strategy:

-   `gh-pages` → live student version
-   `main` → development

When ready: git checkout gh-pages git merge main git push

------------------------------------------------------------------------

# Overall Status

System is stable. Challenge framework works. Phase 2 partially complete.
One full pedagogical challenge implemented. Ready for further refinement
before pushing to production.

Future-me: you're in a good place.

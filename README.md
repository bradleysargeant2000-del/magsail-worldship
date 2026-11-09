# Magnetic sail braking at world-ship scale: code

Code for the paper "Magnetic sail braking at world-ship scale: a size-independent limit, and comparison with fusion braking" (B. Sargeant, 2026).

Requirements: Python 3.10+, NumPy, SciPy. Install with `pip install numpy scipy`.

## Files

| File | What it does |
| --- | --- |
| `full_ship.py` | Habitat mass (Table 1), catalysed D-D fuel physics, rocket-equation mission model with engine and radiator mass, and the trip-time optimiser (Table 2). Inputs are listed at the top (Appendix A of the paper). |
| `contour_model.py` | Unified magnetic sail drag: the field-contour boundary where the loop's own field equals sqrt(mu0 rho v^2), with Fujita (2004) Eq. 12 drag coefficient (3.6 in the fluid limit) and kinetic correction (Section 5). Running it checks the dipole and thin-tube limits and continuity. Set `CDMODEL=ashida` for the 2.9 refit, or `CDSCALE=1.1944` for Fujita's upper bound (4.3). |
| `emergency.py` | Sail as the only brake (Table 4): slowing a 2.82 Mt ship from 4.75% of light speed to 100 km/s, on an 80 x 80 grid of radius and current. |
| `scale_model.py` | `sail_mass()`: superconducting wire, hoop-tension structure and tether mass (Section 6). (Its older `drag()` function is superseded by `contour_model.py`.) |
| `fine_opt.py` | Sail-assisted search on the same fine grid as fusion-only; shows the zero-sail limit reproduces fusion-only results (Section 8). |
| `share_opt.py` | Sail-assisted search where the sail must do at least a fraction f of the braking, for chosen wire and fibre properties (Tables 3 and 5). |
| `sens.py` | Sensitivity runs (Table 6). |
| `recheck_independent.py` | Independent checks: hand-computed habitat, step-by-step flight simulation, hoop tension from inductance, literature reproductions (Freeland/Icarus, Gros). |

## Reproducing the tables

- **Table 1 and Table 2:** `python3 full_ship.py`
- **Tables 3 and 5:** `python3 share_opt.py F J SIGMA RHO NAME`, where F = 0.1 (sail does at least 10% of braking), J = wire current density in A/m^2, SIGMA = allowable fibre stress in Pa (strength / safety factor 2), RHO = fibre density in kg/m^3. Cases used:
  - today: `0.1 3e9 2.9e9 1560`
  - wire 10x: `0.1 3e10 2.9e9 1560`; wire 30x: `0.1 1e11 2.9e9 1560`
  - fibre 3.5x: `0.1 3e9 10e9 1300`
  - wire 5x + fibre: `0.1 1.5e10 10e9 1300`; wire 10x + fibre: `0.1 3e10 10e9 1300`; wire 30x + fibre: `0.1 1e11 10e9 1300`
  - Gas density: set `NDEN=0.03` or `NDEN=0.12` (ions per cm^3; default 0.06).
  - Drag coefficient: default Fujita 3.6; `CDMODEL=ashida` gives 2.9; `CDSCALE=1.1944` gives 4.3.
  - Grid refinement: `RLO=1e7 RHI=1e8 NR=41 ILO=1e6 IHI=5e6 NI=41`.
  Each run takes a few minutes and prints the best design per trip length: departure mass and (radius km, current MA, sail Mt, cruise speed /c, handover speed km/s, sail braking years).
- **Zero-sail check (Section 8):** `python3 fine_opt.py`
- **Table 4 (emergency braking):** `python3 emergency.py`
- **Table 6:** `python3 sens.py`

## Licence

MIT.

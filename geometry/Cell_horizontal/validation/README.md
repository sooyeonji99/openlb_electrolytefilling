# Stability ablations

The production default is preserved: NS BGK, phase BGK with phase tau equal to
liquid NS tau (0.53), legacy internally pinned reservoir. No phi clipping is used.

New diagnostic options:

- `--ns-collision bgk|trt|itrt`: bulk NS collision only. TRT MAGIC=3/16.
  Local TAU_EFF remains the same; density and viscosity inputs do not change.
  The pressure boundary implementation is unchanged.
- `--phase-tau value`: independent conservative Allen-Cahn BGK relaxation time.
  This changes phase mobility; it is not a viscosity-preserving phase-transport
  equivalence. At unchanged dx/dt, tau 0.53→1.0 multiplies mobility by 16.6667.
- `--source-mode legacy|equilibrium|finite|decoupled`: legacy resets pressure and phi moments;
  equilibrium resets the CH equilibrium plus velocity/force history at the liquid
  bath; finite performs the same initial liquid bath setup but no subsequent
  reservoir reset. Finite is a diagnostic with pressure=0, not a continuous supply.
- `--sigma 0`: turn off capillary force in NS as a diagnostic. Phase sharpening and
  contact-angle wall processing remain active; this is not a physical replacement.

`multiphase_diagnostics.csv` additionally records min_rho_lattice and min_tau_eff.
Nonfinite fields, phi bounds, nonpositive mixture density and NS tau<=0.5 stop the
run. The requested speed-amplitude cutoff remains disabled. PASS is a numerical
check to the recorded horizon, not a validation of hydrodynamic accuracy or mass
conservation. Source forcing is nonconservative and a full discrete flux balance
is not implemented.

Run the 5um matrix in fresh output directories:

```bash
python3 validation/run_matrix.py
python3 validation/summarize.py
```

Results: `results/stability_v1/comparison.csv` and `comparison.json`.
Each case records its full configuration. Identical physical circles are used;
material masks are unchanged by the diagnostics/options. All primary comparisons
use dx=5um, dt=0.25us, 3000-step horizon=0.75ms, 1ms linear pressure ramp.
0-pressure cases retain the same liquid source and contact-angle physics unless
explicitly labelled finite or sigma0.

Snapshot pressure verification for finite reservoirs checks finiteness rather
than a prescribed pressure; finite reservoirs do not maintain a pressure bath.
`main_before_stability.cpp` is the unmodified solver before these additions.

Additional installed-example comparisons:

- `itrt`: MultiPhaseIncompressibleInterfaceTRTdynamics, TAU_MINUS=1.5,
  MAXNABLARHO=(rho_lattice_liquid-rho_lattice_gas)/W, as in the installed
  `examples/porousMedia/gasStorage2d` example. No dispersion/phi limiter is added.
- `decoupled`: exclude the prescribed internal bath from Liang force coupling;
  reset bath CH equilibrium/history, NS liquid density/tau and zero gradient/force.
  The surrounding pore coupling, source dimensions and prescribed pressure stay
  unchanged. This tests reservoir consistency, not an external inlet replacement.

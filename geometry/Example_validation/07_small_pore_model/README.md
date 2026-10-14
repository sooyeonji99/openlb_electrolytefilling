# 07 Integrated small-pore model

This final integration check applies the validated conservative phase-field and
incompressible BGK model to a 100 x 100 micrometre domain containing two
12 micrometre throats and two expanded cavities. It uses `dx=2 um`,
`dt=4e-7 s`, `tau_liquid=0.8`, a 60 degree wall angle, and an inlet pressure
1.3 times the nominal throat entry pressure.

After 80,000 steps (`0.032 s` physical time):

- the interface is at `72.284 um`;
- it has crossed the second throat exit at `72.000 um`;
- liquid saturation is `0.41937`;
- maximum lattice speed is `0.001212`, safely below 0.01;
- all recorded fields remain finite.

Status: **PASS**. The run demonstrates stable infiltration through both small
throats and entry into the final cavity. It is a deterministic reduced pore
network and is the direct precursor to importing the 100 x 100 micrometre
random electrode mask.

Run `python3 analyze.py` for the machine-readable result. Open
`results/integrated/tmp/small_pore_model.pvd` in ParaView; the time series and
final saturation are in `small_pore_diagnostics.csv`.

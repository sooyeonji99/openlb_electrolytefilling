"""Check actual MPI histories against the configured on/rest schedule."""
import csv
import math
from pathlib import Path
import sys
from excitation_schedule import active, excitation_at, frequency_at, phase_at

for argument in sys.argv[1:]:
    folder = Path(argument)
    with (folder/'run_config.csv').open() as f:
        config = dict(csv.reader(f))
    with (folder/'multiphase_validation.csv').open() as f:
        final = next(csv.DictReader(f))
    assert final['status'] == 'PASS'
    assert int(final['steps']) == int(config['requested_steps'])
    states = set()
    with (folder/'forcing_history.csv').open() as f:
        for row in csv.DictReader(f):
            time = float(row['time_s'])
            expected = excitation_at(config, time)
            assert math.isclose(float(row['excitation_frequency_hz']), frequency_at(config, time), abs_tol=1e-12)
            assert math.isclose(float(row['phase_rad']), phase_at(config, time), abs_tol=1e-12)
            assert int(row['excitation_active']) == active(config, time)
            assert math.isclose(float(row['inertial_acceleration_y_m_s2']), expected, abs_tol=1e-12)
            assert math.isclose(float(row['effective_acceleration_y_m_s2']), float(config['gravity_y_m_s2']) + expected, abs_tol=1e-12)
            assert float(row['source_pressure_pa']) == 2500
            if time < float(config['excitation_end_s']) - 1e-15:
                states.add(int(row['excitation_active']))
    with (folder/'multiphase_diagnostics.csv').open() as f:
        assert all(row['finite'] == '1' for row in csv.DictReader(f))
    with (folder/'vtk_validation.csv').open() as f:
        rows = list(csv.DictReader(f))
        assert rows and all(row['finite_fluid_fields'] == '1' for row in rows)
    if float(config['excitation_set_period_s']) < float(config['excitation_end_s']):
        assert states == {0, 1}, 'Transition smoke test must exercise both on and rest'
    print(f'{folder.name}: PASS, {final["steps"]} steps, on/rest states={sorted(states)}, pressure/gravity/forcing/VTK checked')

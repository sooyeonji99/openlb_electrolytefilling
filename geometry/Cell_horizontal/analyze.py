"""Validate stored runs and compare matched physical-time pressure experiments."""
import csv,json,math
from pathlib import Path
root=Path(__file__).parent
runs={}
for directory in sorted((root/'results').glob('*')):
    path=directory/'multiphase_diagnostics.csv'
    if not path.exists() or not directory.name.startswith(('verified_', 'mixed_', 'open_', 'short_')):continue
    rows=list(csv.DictReader(path.open()))
    if not rows:continue
    finite=all(all(math.isfinite(float(row[k])) for k in ('time_s','pore_mean_phi','min_phi','max_phi','max_speed_lattice')) and row['finite']=='1' for row in rows)
    config=dict(csv.reader((directory/"run_config.csv").open()))
    first,last=rows[0],rows[-1]
    runs[directory.name]={'outer_boundary':config.get('outer_boundary','closed'),'geometry_version':config.get('geometry_version','legacy_single_diameter_v1'),'checks_pass':finite,'time_s':float(last['time_s']),'steps':int(last['step']),'initial_saturation':float(first['pore_mean_phi']),'final_saturation':float(last['pore_mean_phi']),'max_speed_lattice':max(float(r['max_speed_lattice']) for r in rows)}
for directory in sorted((root/'results').iterdir()):
    if not directory.name.startswith(('closed_cell_', 'mixed_')):continue
    summary=directory/'closed_cell_analysis.json'
    if summary.exists():runs[directory.name]=json.loads(summary.read_text())
validation_summary=root/'results/stability_v1/comparison.json'
if validation_summary.exists():runs['stability_validation']={'report':'results/stability_v1/REPORT.md','cases':json.loads(validation_summary.read_text())}
(root/'results'/'analysis.json').write_text(json.dumps(runs,indent=2)+'\n')
print(json.dumps(runs,indent=2))

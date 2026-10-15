"""Collect all completed and running stability ablations without altering results."""
import csv,json,math
from pathlib import Path
root=Path(__file__).resolve().parent.parent
base=root/'results/stability_v1';summary=[]
for folder in sorted(base.iterdir()):
 if not folder.is_dir():continue
 path=folder/'multiphase_diagnostics.csv'
 if not path.exists():continue
 rows=list(csv.DictReader(path.open()))
 if not rows:continue
 config=dict(csv.reader((folder/'run_config.csv').open()));last=rows[-1]
 result_path=folder/'multiphase_validation.csv';result=next(csv.DictReader(result_path.open())) if result_path.exists() else {}
 entry={'case':folder.name,'status':result.get('status','RUNNING'),'dx_um':float(config['dx_m'])*1e6,'steps':int(last['step']),'time_ms':float(last['time_s'])*1000,'applied_pressure_pa':float(last['reservoir_pressure_pa']),'min_phi':float(last['min_phi']),'max_phi':float(last['max_phi']),'min_rho_lattice':float(last['min_rho_lattice']),'min_tau_eff':float(last['min_tau_eff']),'max_speed_over_run':max(float(r['max_speed_lattice']) for r in rows),'ns_collision':config['ns_collision'],'source_mode':config['source_mode'],'phase_tau':float(config['phase_tau']),'sigma':float(config['sigma_N_m'])}
 entry['min_phi_over_run']=min(float(r['min_phi']) for r in rows)
 entry['max_phi_over_run']=max(float(r['max_phi']) for r in rows)
 for threshold in (.1,.3,1):
  first=next((r for r in rows if float(r['max_speed_lattice'])>=threshold),None);entry['first_speed_'+str(threshold)+'_step']=int(first['step']) if first else None
 summary.append(entry)
(base/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
with (base/'comparison.csv').open('w') as f:
 if summary:w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
print(json.dumps(summary,indent=2))

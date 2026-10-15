"""One-factor stability ablations; separate output directories; no phi clipping."""
import concurrent.futures,csv,json,subprocess,time
from pathlib import Path
root=Path(__file__).resolve().parent.parent
base=root/'results/stability_v1';base.mkdir(exist_ok=True)
if any((base/name).exists() for name in ['bgk_p0','bgk_ramp','trt_ramp']):raise SystemExit('Existing matrix results; choose a new base directory.')
cases=[
 ('bgk_p0',{'pressure':0}),
 ('bgk_ramp',{}),
 ('trt_ramp',{'collision':'trt'}),
 ('source_equilibrium_ramp',{'source':'equilibrium'}),
 ('finite_p0',{'source':'finite','pressure':0}),
 ('phase_tau1_ramp',{'phase':1}),
 ('bgk_p0_sigma0',{'pressure':0,'sigma':0}),
]
def run(item):
 name,o=item;folder=base/name
 cmd=[str(root/'main'),'--dx-um','5','--steps','3000','--diagnostic-every','1','--output-every','1000','--reservoir-pressure-pa',str(o.get('pressure',10000)),'--pressure-ramp-s','.001','--ns-collision',o.get('collision','bgk'),'--source-mode',o.get('source','legacy'),'--phase-tau',str(o.get('phase',.53)),'--sigma',str(o.get('sigma',1e-5)),'--output',str(folder)]
 start=time.monotonic()
 with (base/(name+'.log')).open('w') as log:result=subprocess.run(cmd,cwd=root,stdout=log,stderr=subprocess.STDOUT)
 summary={'case':name,'returncode':result.returncode,'wall_seconds':time.monotonic()-start,'command':cmd}
 path=folder/'multiphase_validation.csv'
 if path.exists():summary.update(next(csv.DictReader(path.open())))
 rows=list(csv.DictReader((folder/'multiphase_diagnostics.csv').open())) if (folder/'multiphase_diagnostics.csv').exists() else []
 if rows:summary.update(final_diagnostics=rows[-1],max_speed_over_run=max(float(r['max_speed_lattice']) for r in rows))
 (folder/'experiment_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 print(name,summary.get('status'),summary.get('steps'),flush=True)
 return summary
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 result=list(pool.map(run,cases))
(base/'matrix_summary.json').write_text(json.dumps(result,indent=2)+'\n')

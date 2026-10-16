"""Two isolated 0 Pa controls on identical geometry; no production changes."""
import concurrent.futures,csv,hashlib,json,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
BASE=ROOT/'results/zero_pressure_controls_v1'

def run(case):
    name,sigma=case;folder=BASE/name
    if folder.exists():raise RuntimeError(f'Refusing to overwrite {folder}')
    cmd=[str(ROOT/'main'),'--dx-um','5','--tau','.53','--phase-tau','.53',
         '--nu-liquid','1e-6','--nu-gas','1e-7','--outer-boundary','open',
         '--ns-collision','bgk','--source-mode','legacy',
         '--reservoir-pressure-pa','0','--pressure-ramp-s','.00025','--sigma',str(sigma),
         '--steps','10000','--diagnostic-every','1','--output-every','100','--output',str(folder)]
    start=time.monotonic()
    with (BASE/(name+'.log')).open('w') as stream:
        proc=subprocess.run(cmd,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
    summary={'case':name,'sigma_N_m':sigma,'command':cmd,'returncode':proc.returncode,'solver_wall_seconds':time.monotonic()-start}
    validation=folder/'multiphase_validation.csv'
    if validation.exists():
        with validation.open() as stream:summary.update(next(csv.DictReader(stream)))
    diagnostics=folder/'multiphase_diagnostics.csv'
    if diagnostics.exists():
        with diagnostics.open() as stream:rows=list(csv.DictReader(stream))
        summary.update(last_diagnostic=rows[-1],max_speed_over_run=max(float(r['max_speed_lattice']) for r in rows),min_phi_over_run=min(float(r['min_phi']) for r in rows),max_phi_over_run=max(float(r['max_phi']) for r in rows),diagnostic_rows=len(rows),diagnostic_continuous=all(int(r['step'])==i for i,r in enumerate(rows)))
    with (BASE/(name+'_vtk_check.log')).open('w') as stream:
        check=subprocess.run(['python3',str(ROOT/'verify_vtk.py'),str(folder)],cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
    summary['vtk_check_returncode']=check.returncode
    summary['circles_sha256']=hashlib.sha256((folder/'circles.csv').read_bytes()).hexdigest()
    summary['full_10000_pass']=proc.returncode==0 and summary.get('status')=='PASS' and int(summary.get('steps',0))==10000 and check.returncode==0 and summary.get('diagnostic_continuous',False)
    (folder/'experiment_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(name,summary.get('status'),summary.get('steps'),'maxU',summary.get('max_speed_over_run'),flush=True)
    return summary

if __name__=='__main__':
    BASE.mkdir(parents=True,exist_ok=True)
    if any((BASE/name).exists() for name in ('p0_sigma_default','p0_sigma0')):raise SystemExit('Existing control output; select a new BASE')
    (BASE/'source_sha256.json').write_text(json.dumps({name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['main.cpp','geometry.h','run_controls.py']},indent=2)+'\n')
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        summaries=[]
        futures=[pool.submit(run,x) for x in [('p0_sigma_default',1e-5),('p0_sigma0',0)]]
        for future in concurrent.futures.as_completed(futures):
            summaries.append(future.result());(BASE/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    (BASE/'completion.json').write_text(json.dumps({'complete':True,'same_circles':len({s['circles_sha256'] for s in summaries})==1,'cases':len(summaries)},indent=2)+'\n')
    print('ALL CONTROLS COMPLETE',flush=True)

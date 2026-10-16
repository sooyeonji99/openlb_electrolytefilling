"""2x2 controls: contact angle (30/90 deg) and gas nu (1e-7/1e-6).
Fixed: geometry, rho, sigma, liquid viscosity, dt, NS liquid and phase tau.
The gas change simultaneously changes its physical viscosity and NS tau;
it is a sensitivity control, not an isolated collision-algorithm comparison.
"""
import concurrent.futures,csv,hashlib,json,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
BASE=ROOT/'results/wetting_gas_theta30_v1'
CASES=[('theta30_gas_baseline',30,1e-7),('theta90_gas_baseline',90,1e-7),('theta30_gas_matched',30,1e-6),('theta90_gas_matched',90,1e-6)]

def run(case):
    name,theta,nu_gas=case;folder=BASE/name
    if folder.exists():raise RuntimeError(f'Refusing to overwrite {folder}')
    cmd=[str(ROOT/'main'),'--contact-angle-deg',str(theta),'--dx-um','5','--tau','.53','--phase-tau','.53',
         '--phase-projection','off','--material-interpolation','linear','--rho-gas','100','--nu-liquid','1e-6','--nu-gas',str(nu_gas),'--outer-boundary','open',
         '--ns-collision','bgk','--source-mode','legacy','--sigma','1e-5',
         '--reservoir-pressure-pa','0','--pressure-ramp-s','.00025',
         '--steps','10000','--diagnostic-every','1','--output-every','100','--output',str(folder)]
    start=time.monotonic()
    with (BASE/(name+'.log')).open('w') as f:proc=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    r={'case':name,'contact_angle_deg':theta,'nu_gas_m2_s':nu_gas,'command':cmd,'returncode':proc.returncode,'solver_wall_seconds':time.monotonic()-start}
    with (folder/'multiphase_validation.csv').open() as f:r.update(next(csv.DictReader(f)))
    with (folder/'multiphase_diagnostics.csv').open() as f:rows=list(csv.DictReader(f))
    r.update(last_diagnostic=rows[-1],max_speed_over_run=max(float(x['max_speed_lattice']) for x in rows),min_phi_over_run=min(float(x['min_phi']) for x in rows),max_phi_over_run=max(float(x['max_phi']) for x in rows),min_rho_over_run=min(float(x['min_rho_lattice']) for x in rows),min_ns_tau_over_run=min(float(x['min_tau_eff']) for x in rows),diagnostic_rows=len(rows),diagnostic_continuous=all(int(x['step'])==i for i,x in enumerate(rows)))
    r['first_speed_above_0_1']=next((int(x['step']) for x in rows if float(x['max_speed_lattice'])>.1),None)
    with (BASE/(name+'_vtk_check.log')).open('w') as f:check=subprocess.run(['python3',str(ROOT/'verify_vtk.py'),str(folder)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    r['vtk_check_returncode']=check.returncode
    r['circles_sha256']=hashlib.sha256((folder/'circles.csv').read_bytes()).hexdigest()
    r['full_10000_pass']=proc.returncode==0 and r['status']=='PASS' and int(r['steps'])==10000 and check.returncode==0 and r['diagnostic_continuous'] and len(rows)==10001
    (folder/'experiment_summary.json').write_text(json.dumps(r,indent=2)+'\n')
    print(name,r['status'],r['steps'],'maxU',r['max_speed_over_run'],flush=True)
    if proc.returncode not in (0,2):raise RuntimeError(f'Execution failed: {name}')
    return r

if __name__=='__main__':
    BASE.mkdir(parents=True,exist_ok=True)
    if any((BASE/c[0]).exists() for c in CASES):raise SystemExit('Existing results; select a new BASE')
    snapshot=BASE/'source_at_execution';snapshot.mkdir()
    hashes={}
    for name in ['main.cpp','geometry.h','run_wetting_gas_checks.py','verify_vtk.py']:
        data=(ROOT/name).read_bytes();(snapshot/name).write_bytes(data);hashes[name]=hashlib.sha256(data).hexdigest()
    (BASE/'source_sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
    plan={'cases':CASES,'steps':10000,'pressure_pa':0,'sigma_N_m':1e-5,'phase_tau':.53,'liquid_NS_tau':.53,'gas_NS_tau_baseline':.503,'gas_NS_tau_matched':.53,'note':'90 deg is a neutral wetting control; wall boundary remains active. Matched gas viscosity is an artificial sensitivity case; rho_gas stays 100 kg/m3.'}
    (BASE/'PLAN.json').write_text(json.dumps(plan,indent=2)+'\n')
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(run,case) for case in CASES];summaries=[]
        for future in concurrent.futures.as_completed(futures):
            summaries.append(future.result());(BASE/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    (BASE/'completion.json').write_text(json.dumps({'complete':True,'cases':len(summaries),'same_circles':len({r['circles_sha256'] for r in summaries})==1},indent=2)+'\n')
    print('ALL WETTING/GAS CHECKS COMPLETE',flush=True)

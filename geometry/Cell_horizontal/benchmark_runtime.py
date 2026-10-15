"""Measure current solver throughput including per-step diagnostics, excluding startup/VTK."""
import json,subprocess,time
from pathlib import Path
root=Path(__file__).resolve().parent
folder=root/'results/short_dx02_benchmark100'
if folder.exists():folder=folder.with_name(folder.name+'_'+str(time.time_ns()))
log=folder.with_suffix('.log')
command=[str(root/'main'),'--steps','100','--diagnostic-every','1','--output-every','100','--output',str(folder)]
samples=[];last=-1;start=time.perf_counter()
with log.open('w') as output:
    process=subprocess.Popen(command,cwd=root,stdout=output,stderr=subprocess.STDOUT)
    while True:
        path=folder/'multiphase_diagnostics.csv'
        if path.exists():
            lines=path.read_text().splitlines()
            if len(lines)>1:
                try:step=int(lines[-1].split(',')[0])
                except ValueError:step=-1
                if step>last:samples.append((step,time.perf_counter()-start));last=step
        if process.poll() is not None:break
        time.sleep(.05)
wall=time.perf_counter()-start
middle=[(step,t) for step,t in samples if 10<=step<=90]
if process.returncode!=0 or len(middle)<2:raise SystemExit('Benchmark incomplete; inspect '+str(log))
first,final=middle[0],middle[-1]
seconds_per_step=(final[1]-first[1])/(final[0]-first[0])
result={'command':command,'returncode':process.returncode,'benchmark_steps':100,'total_wall_seconds_including_initialization_and_vtk':wall,'timed_first_step':first[0],'timed_last_step':final[0],'seconds_per_step':seconds_per_step,'steps_per_second':1/seconds_per_step,'estimated_30000_compute_seconds':30000*seconds_per_step,'estimated_30000_compute_hours':30000*seconds_per_step/3600,'conditions':'one solver process, default dx=2um, per-step diagnostics; output cost excluded from throughput; not a stability prediction'}
(folder/'runtime_benchmark.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

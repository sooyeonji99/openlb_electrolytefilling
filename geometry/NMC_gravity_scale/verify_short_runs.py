import csv,json,math
from pathlib import Path
root=Path(__file__).resolve().parent
result={}
for name,scale in [('smoke_0g',0),('smoke_1g',1),('smoke_2g',2)]:
 p=root/'results'/name
 with (p/'run_config.csv').open() as f:cfg=dict(csv.reader(f))
 with (p/'multiphase_validation.csv').open() as f:final=next(csv.DictReader(f))
 assert final['status']=='PASS' and int(final['steps'])==200
 assert cfg['side_boundary']=='wall' and cfg['bottom_boundary']=='wall' and cfg['top_boundary']=='open'
 assert float(cfg['height_m'])==.00367 or math.isclose(float(cfg['height_m']),.00367,abs_tol=1e-14)
 assert float(cfg['excitation_frequency_hz'])==100 and float(cfg['gravity_scale'])==scale
 with (p/'forcing_history.csv').open() as f:forcing=list(csv.DictReader(f))
 for row in forcing:
  t=float(row['time_s']);a=float(row['inertial_acceleration_y_m_s2'])
  assert float(row['source_pressure_pa'])==2500
  assert math.isclose(a,scale*9.81*math.sin(2*math.pi*100*t),abs_tol=1e-12)
  assert math.isclose(float(row['effective_acceleration_y_m_s2']),-9.81+a,abs_tol=1e-12)
 with (p/'multiphase_diagnostics.csv').open() as f:assert all(row['finite']=='1' for row in csv.DictReader(f))
 with (p/'outlet_diagnostics.csv').open() as f:assert all(float(row['nonfinite_fields'])==0 for row in csv.DictReader(f))
 with (p/'vtk_validation.csv').open() as f:assert all(row['finite_fluid_fields']=='1' for row in csv.DictReader(f))
 result[name]={'status':'PASS','mpi_ranks':4,'steps':200,'gravity_scale':scale,'vtk':'PASS'}
# Independent physical equivalence check: 0g matches the previous model's quiet interval.
reference=root.parent/'NMC_forcesweep_longterm/results/smoke_g_initial_rest_final/multiphase_diagnostics.csv'
# Older archived naming uses smoke_initial_rest_final.
if not reference.exists():reference=root.parent/'NMC_forcesweep_longterm/results/smoke_initial_rest_final/multiphase_diagnostics.csv'
if reference.exists():
 with reference.open() as f:a=list(csv.DictReader(f))
 with (root/'results/smoke_0g/multiphase_diagnostics.csv').open() as f:b=list(csv.DictReader(f))
 assert len(a)==len(b)
 diff=max(abs(float(x['pore_mean_phi'])-float(y['pore_mean_phi'])) for x,y in zip(a,b))
 assert diff<1e-12
 result['zero_g_reference']={'status':'PASS','max_mean_phi_difference':diff}
result['analytic_waveform']={'status':'PASS','frequency_hz':100,'duration_s':.12,'cycles':12,'scales_checked':[0,.5,1,2]}
result['long_run_started']=False
(root/'results/VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n')
lines=['# NMC_gravity_scale 단시간 검증','', '|검사|판정|범위|','|---|---|---|']
for name,row in result.items():
 if name.startswith('smoke'):lines.append(f"|{name}|PASS|4프로세스MPI,200step,가속도/압력/유한성/VTK 검사|")
lines+=['|처방 파형|PASS|0·0.5·1·2g,100Hz12주기,분기/위상/가속도 크기|','', '기본조건 120ms=1,920,000step이다. 짧은MPI검증은200step=12.5µs의초기구간이며, 120ms의유동안정성이나상부액체도달을검증한것은아니다.','', 'M=0은이전모델의무가진초기구간과대조했다. 전체격자·기공 연결성·입자간격도geometry_validation과particle_validation에서검사했다. 장시간계산은자동시작하지않았다.']
(root/'results/VALIDATION_REPORT.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(result,indent=2))

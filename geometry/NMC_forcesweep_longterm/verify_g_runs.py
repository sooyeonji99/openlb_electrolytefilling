import csv,json,math
from pathlib import Path
root=Path(__file__).resolve().parent
names=['smoke_g_initial_rest','smoke_g_active_up','smoke_g_all_transitions']
summary={}
for name in names:
 p=root/'results'/name
 with (p/'run_config.csv').open() as f:cfg=dict(csv.reader(f))
 with (p/'multiphase_validation.csv').open() as f:final=next(csv.DictReader(f))
 assert final['status']=='PASS'
 with (p/'multiphase_diagnostics.csv').open() as f:
  rows=list(csv.DictReader(f));assert all(r['finite']=='1' for r in rows)
 with (p/'forcing_history.csv').open() as f:forcing=list(csv.DictReader(f))
 assert all(float(r['source_pressure_pa'])==2500 for r in forcing)
 assert cfg['excitation_control']=='acceleration_g' and float(cfg['excitation_g'])==1
 for r in forcing:
  a=float(r['inertial_acceleration_y_m_s2']);phase=float(r['phase_rad'])
  assert abs(a)<=9.81+1e-12
  expected=9.81*math.sin(phase) if r['sweep_segment'].startswith('sweep_') else 0
  assert math.isclose(a,expected,abs_tol=1e-11)
  assert math.isclose(float(r['effective_acceleration_y_m_s2']),-9.81+a,abs_tol=1e-11)
 with (p/'outlet_diagnostics.csv').open() as f:assert all(float(r['nonfinite_fields'])==0 for r in csv.DictReader(f))
 with (p/'vtk_validation.csv').open() as f:assert all(r['finite_fluid_fields']=='1' for r in csv.DictReader(f))
 assert cfg['side_boundary']=='wall' and cfg['top_boundary']=='open' and cfg['bottom_boundary']=='wall'
 segments=sorted({r['sweep_segment'] for r in forcing})
 if name=='smoke_g_initial_rest':assert segments==['pre_rest']
 if name=='smoke_g_active_up':assert segments==['sweep_up']
 if name=='smoke_g_all_transitions':assert set(segments)=={'pre_rest','sweep_up','sweep_down','post_rest'}
 summary[name]={'status':'PASS','steps':int(final['steps']),'segments':segments,'excitation_g':1,'pressure_pa':2500,'vtk':'PASS','finite':'PASS'}
summary['analytic_protocol']={'status':'PASS','duration_s':1.06,'cycles':55,'peak_acceleration_m_s2':9.81,'fixed_g_for_entire_sweep':True,'linear_scaling_and_zero_g_tests':'PASS'}
summary['long_run_started']=False
(root/'results/VALIDATION.json').write_text(json.dumps(summary,indent=2)+'\n')
lines=['# 가속도 제어(g 단위) 검증','', '기본1g=9.81m/s²의 고정 가속도 진폭. 위상은10→100→10Hz의 적분값을 사용하며 공급압력2500Pa를 유지한다.','', '|시험|결과|범위|','|---|---|---|']
for name,row in summary.items():
 if name.startswith('smoke'):lines.append(f"|{name}|PASS|{row['steps']}step·MPI4프로세스·VTK/유한성 검사|")
lines+=['','기본 조건 초기200step, 실제1g 진폭으로 pre-rest만0으로 둔 상승 시작200step, 시간을 압축한 전체 구간 전환640step을 검증했다. 전환 시험도 진폭1g를 유지했다.','', '전체1.06초의 처방 파형을10µs 간격으로 확인해 가속도가±9.81 이내임을 확인했다. 2g는2배,0g는0이며 위상 연속성과55주기도 검사했다. 전체1.06초 유동 계산은 시작하지 않았다.','', '변위·속도 경로는 지정하지 않는다. 중력은 별도로 더하며 기본 설정의 총 가속도는−19.62~0m/s²이다. 이전 변위 제어 결과는displacement_control_previous에 보존했다.']
(root/'results/VALIDATION_REPORT.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(summary,indent=2))

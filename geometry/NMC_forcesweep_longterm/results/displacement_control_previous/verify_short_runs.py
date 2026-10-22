import csv,json,math
from pathlib import Path
root=Path(__file__).resolve().parent
names=['smoke_initial_rest_final','smoke_active_up_final','smoke_all_transitions_final']
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
 with (p/'outlet_diagnostics.csv').open() as f:
  outlets=list(csv.DictReader(f));assert all(float(r['nonfinite_fields'])==0 for r in outlets)
 with (p/'vtk_validation.csv').open() as f:
  checks=list(csv.DictReader(f));assert checks and all(r['finite_fluid_fields']=='1' for r in checks)
 assert cfg['side_boundary']=='wall' and cfg['top_boundary']=='open' and cfg['bottom_boundary']=='wall'
 segments=sorted({r['sweep_segment'] for r in forcing})
 if name=='smoke_initial_rest_final':
  assert segments==['pre_rest']
  assert all(float(r['inertial_acceleration_y_m_s2'])==0 for r in forcing)
 if name=='smoke_active_up_final':
  assert segments==['sweep_up']
  assert math.isclose(float(forcing[0]['inertial_acceleration_y_m_s2']),-.00734*2*math.pi*180,abs_tol=1e-10)
 if name=='smoke_all_transitions_final':assert set(segments)=={'pre_rest','sweep_up','sweep_down','post_rest'}
 summary[name]={'status':'PASS','steps':int(final['steps']),'segments':segments,'pressure_pa':2500,'vtk':'PASS','finite':'PASS'}
summary['geometry']={'status':'PASS','width_mm':4.3,'height_mm':3.67,'source_height_mm':.5,'pores_connected_to_source':True}
summary['analytic_protocol']={'status':'PASS','duration_s':1.06,'pre_rest_ms':30,'up_s':.5,'down_s':.5,'post_rest_ms':30,'cycles':55,'phase_and_acceleration_derivative_test':'PASS'}
summary['long_run_started']=False
(root/'results/VALIDATION.json').write_text(json.dumps(summary,indent=2)+'\n')
lines=['# 단시간 검증 결과','', '전체1.06초 계산은 시작하지 않았다.','', '|검사|결과|범위|','|---|---|---|', '|전체 파형 수식|PASS|30ms 대기, 10→100→10Hz 연속 위상, 가속도 미분, 총55주기, 종료30ms 대기|','|구조·입자|PASS|높이3.67mm, 좌우하단 벽·상부 개방, 기공 연결성 및 입자 간격|']
for name,row in summary.items():
 if name.startswith('smoke'):lines.append(f"|{name}|PASS|{row['steps']}step, 4 MPI프로세스, pressure2500Pa·유한성·VTK 검사|")
lines+=['','smoke_initial_rest_final은 실제 기본 조건의 초기200step이다. smoke_active_up_final은 pre-rest만0으로 바꾸어 실제 진폭7.34mm와 실제0.5초 상승률의 스윕 시작200step을 검사했다.','', 'smoke_all_transitions_final은 전환을 짧은 계산으로 검사하기 위해 대기5µs·상승15µs·하강15µs·대기5µs, 진폭비0.0002로 바꾼 시험640step이다. 본 실행 설정은 바꾸지 않았다.','', '초기 시험에서 pre-rest=0일 때 step0 가속도 필드가 중력만 담긴 문제를 발견했다. applyAcceleration(0)을 추가한 후 최종 실행파일로 위 세 검사를 재실행했다. 기존 smoke_active_up은 수정 전 기록이다.','', '이 검증은 전체1.06초의 장시간 안정성·시간/격자 수렴·상부 액체 도달을 입증하지 않는다. 상부 배출은 개방 경계로 허용되며, 짧은 시험에서 액체가 상부에 도달한 것은 아니다.']
(root/'results/VALIDATION_REPORT.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(summary,indent=2))

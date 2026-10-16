"""Generate the report after both 10000-step zero-pressure controls finish."""
import csv,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
BASE=ROOT/'results/zero_pressure_controls_v1'
complete=json.loads((BASE/'completion.json').read_text())
assert complete['complete'] and complete['cases']==2 and complete['same_circles']
rows=json.loads((BASE/'summary.json').read_text());cases={r['case']:r for r in rows}
original=ROOT/'results/geometry_dx05/circles.csv'
checks={}
for name,r in cases.items():
    folder=BASE/name
    cfg=dict(csv.reader((folder/'run_config.csv').open()))
    assert float(cfg['pressure_pa'])==0 and float(cfg['tau_liquid'])==.53 and float(cfg['phase_tau'])==.53
    assert float(cfg['sigma_N_m'])==r['sigma_N_m']
    assert r['diagnostic_continuous'] and r['diagnostic_rows']==int(r['steps'])+1
    assert (folder/'circles.csv').read_bytes()==original.read_bytes()
    for n in ('geometry_validation.csv','particle_validation.csv'):
        with (folder/n).open() as f:assert next(csv.DictReader(f))['status']=='PASS'
    with (folder/'vtk_validation.csv').open() as f:vtks=list(csv.DictReader(f))
    assert r['vtk_check_returncode']==0 and all(x['finite_fluid_fields']=='1' for x in vtks)
    checks[name]={'geometry_identical':True,'saved_vtk_checked':len(vtks),'input_parameters_confirmed':True}
hashes=json.loads((BASE/'source_sha256.json').read_text())
snapshot=BASE/'source_at_execution'
assert all(hashlib.sha256(((snapshot if snapshot.exists() else ROOT)/name).read_bytes()).hexdigest()==digest for name,digest in hashes.items())
(BASE/'independent_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
keys=['case','status','steps','physical_time_s','sigma_N_m','max_speed_over_run','min_phi_over_run','max_phi_over_run','solver_wall_seconds','vtk_check_returncode']
with (BASE/'comparison.csv').open('w') as f:
    w=csv.DictWriter(f,keys,extrasaction='ignore');w.writeheader();w.writerows(rows)
lines=['# 0 Pa 대조군 및 고체 모델 확인','',
       '두 실행 모두 동일한 확대 기공 형상, 5 µm 격자, NS/상장 tau=0.53, 기본 액체/기체 물성, legacy 공급원과 개방 외곽 조건입니다.',
       '외부 압력 차는 0 Pa로 고정했습니다. 각 실행은 독립 초기화했으며 원형 입자 CSV는 원본 및 서로 간에 완전히 같습니다.',
       '', '| 조건 | 결과 | 스텝 | 물리 시간 (ms) | 전체 구간 최대 격자 속도 (NS 통계) | 최종 phi 범위 |',
       '|---|---|---:|---:|---:|---|']
for name in ['p0_sigma_default','p0_sigma0']:
    r=cases[name];last=r['last_diagnostic']
    lines.append(f"| {name}, sigma={r['sigma_N_m']:g} N/m | {r['status']} | {r['steps']} | {float(r['physical_time_s'])*1000:g} | {r['max_speed_over_run']:.6g} | {float(last['min_phi']):.6g} ~ {float(last['max_phi']):.6g} |")
lines += ['', '## 해석','',
          '표면장력을 유지한 0 Pa 조건은 1752스텝(0.438 ms)에 phi 하한 이탈, 음의 혼합 밀도 및 NS tau<=0.5로 중단됐습니다.',
          '최종 min_phi=-0.130723, min_rho=-0.0176504, min_NS_tau=0.499470. 속도는 그 전에 증가했습니다.',
          '최종 공간장의 최소 phi는 약 (2.680, 0.205) mm, 최대 속도는 약 (2.705, 0.230) mm의 내부 기공에서 관찰됐습니다.',
          '고압 공급이 없어도 실패하므로 압력 구동은 불안정 발생의 필수 조건이 아닙니다. 2500 Pa 실행과는 계면 경로/상태가 달라 실패 시점만으로 압력 안정화 효과를 일반화하지 않습니다.',
          '',
          '표면장력 0 대조군은 상장 sharpening, wetting 경계, 매 스텝 reservoir phi=1 재설정을 유지합니다.',
          '제거한 것은 NS 계면력에 해당하며 완전히 모든 계면 물리를 제거한 대조군이 아닙니다.',
          '따라서 이 대조군의 안정성은 계면력-유동 피드백의 관여를 지지하지만 계면력 공식의 구현 오류나 특정 벽/충돌 모델 하나를 확정하지 않습니다.',
          '상장 mobility/격자 해상도 및 기체 NS tau=0.503, 고체 wetting, 강제 공급원 처리는 여전히 후속 원인 분리 대상입니다.',
          '', '## 고체 입력값','',
          '고정·비변형의 2차원 원형 벽, 정지 no-slip/불투과 NS bounce-back, 곡면 상장 wetting 경계, 명목 접촉각 60°입니다.',
          '실제 입자 330개, 직경 50.1651–199.310 µm, 최소 간격 약 20 µm, 공급원 제외 공극률 0.398094입니다.',
          '고체 밀도, 탄성계수, 포아송비, 열물성, 전기화학 물성은 설정하지 않았고 고체 역학/열 방정식도 없습니다.',
          '1000/100 kg/m³는 액체/기체 밀도이고 sigma=1e-5 N/m는 액체–기체 표면장력입니다. 실제 NMC 물성을 넣은 모델이 아닙니다.',
          '전체 목록: ../../SOLID_PROPERTIES.md (폴더 루트의 문서).',
          '', '## 판정 범위와 다음 순서','',
          'NS 계면력이 없는 조건의 PASS를 실제 침투 해석 조건으로 채택하지 않습니다. 이번 비교는 원인 분리를 위한 진단입니다.',
          '다음 단계는 같은 0 Pa/기본 sigma에서 중립 젖음(90°) 및 기체 충돌/점성 대조군으로 고체 표면 효과와 기체 유동 감쇠를 분리하는 것입니다.',
          '실제 접촉각/기체 점성의 수정은 물리 조건 변경이므로 안정성 조치와 물성 검증을 구분해야 합니다.',
          '질량 수지, 접촉각 재현, Laplace 압력 및 정적 계면의 기생 유동을 검증한 후 압력 공급 침투로 복귀하는 것이 좋습니다.',
          'main.cpp/geometry.h 및 생산 기본값은 변경하지 않았습니다. 코드 해시와 VTK/입력 일치 검사는 별도 JSON에 기록했습니다.',
          '', '분석 재생성: python3 analyze_controls.py. 실제 실행 명령은 각 experiment_summary.json에 있습니다.']
(BASE/'REPORT.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(checks,indent=2))

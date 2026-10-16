"""Validate the completed 2x2 matrix and generate a bounded-scope report."""
import csv,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
BASE=ROOT/'results/wetting_gas_theta30_v1'
complete=json.loads((BASE/'completion.json').read_text())
assert complete['complete'] and complete['cases']==4 and complete['same_circles']
rows=json.loads((BASE/'summary.json').read_text());cases={r['case']:r for r in rows}
hashes=json.loads((BASE/'source_sha256.json').read_text())
for name,digest in hashes.items():assert hashlib.sha256((BASE/'source_at_execution'/name).read_bytes()).hexdigest()==digest
reference=ROOT/'results/geometry_dx05/circles.csv';checks={}
for name,r in cases.items():
    folder=BASE/name;cfg=dict(csv.reader((folder/'run_config.csv').open()))
    assert float(cfg['contact_angle_deg'])==r['contact_angle_deg']
    assert float(cfg['nu_gas_m2_s'])==r['nu_gas_m2_s']
    assert abs(float(cfg['tau_gas'])-(.5+.03*r['nu_gas_m2_s']/1e-6))<1e-12
    for key,value in [('phase_tau',.53),('tau_liquid',.53),('pressure_pa',0),('sigma_N_m',1e-5),('nu_liquid_m2_s',1e-6),('rho_gas_kg_m3',100),('dx_m',5e-6),('dt_s',2.5e-7)]:assert abs(float(cfg[key])-value)<max(abs(value)*1e-12,1e-18),(name,key)
    assert (folder/'circles.csv').read_bytes()==reference.read_bytes()
    assert r['diagnostic_continuous'] and r['diagnostic_rows']==int(r['steps'])+1
    for filename in ['geometry_validation.csv','particle_validation.csv']:
        with (folder/filename).open() as f:assert next(csv.DictReader(f))['status']=='PASS'
    with (folder/'vtk_validation.csv').open() as f:vtks=list(csv.DictReader(f))
    assert r['vtk_check_returncode']==0 and all(v['finite_fluid_fields']=='1' for v in vtks)
    checks[name]={'geometry_identical':True,'parameters_checked':True,'vtk_snapshots_checked':len(vtks),'diagnostic_rows':r['diagnostic_rows']}
(BASE/'independent_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
keys=['case','contact_angle_deg','nu_gas_m2_s','status','steps','physical_time_s','max_speed_over_run','min_phi_over_run','max_phi_over_run','min_rho_over_run','min_ns_tau_over_run','first_speed_above_0_1','full_10000_pass','solver_wall_seconds','vtk_check_returncode']
with (BASE/'comparison.csv').open('w') as f:
    writer=csv.DictWriter(f,keys,extrasaction='ignore');writer.writeheader();writer.writerows(sorted(rows,key=lambda r:(r['nu_gas_m2_s'],r['contact_angle_deg'])))
lines=['# 접촉각 30° / 젖음·기체 점성 대조군','',
       '현재 생산 기본 접촉각을 60°에서 30°로 변경했습니다. 새 CLI는 --contact-angle-deg입니다.',
       '기존 60° 출력은 변경하지 않았고 실행 당시 소스도 해당 결과 폴더에 해시 검증 후 보관했습니다.',
       '', '## 고정 조건','',
       '확대 기공 형상(공극률 0.398094, 최소 입자 간격 약 20 µm), dx=5 µm, 0 Pa, sigma=1e-5 N/m, NS BGK, legacy 공급원, 하단 벽/상단·좌우 개방.',
       '액체 동점성 1e-6 m²/s, 밀도 1000/100 kg/m³, 액체 NS tau=상장 tau=0.53, dt=0.25 µs, 계면 두께 4격자는 모든 조건에서 같습니다.',
       '각 독립 초기화 실행의 목표는 10000스텝=2.5 ms입니다. 매 스텝 진단과 100스텝 공간 출력을 저장했습니다.',
       '', '## 결과','',
       '| 접촉각 | 기체 동점성 (m²/s) | 기체 NS tau | 결과 | 스텝 | 시간 (ms) | 최대 격자 속도 (NS 통계) | 최종 phi 범위 |',
       '|---:|---:|---:|---|---:|---:|---:|---|']
for r in sorted(rows,key=lambda r:(r['nu_gas_m2_s'],r['contact_angle_deg'])):
    d=r['last_diagnostic'];tau=.5+.03*r['nu_gas_m2_s']/1e-6
    lines.append(f"| {r['contact_angle_deg']}° | {r['nu_gas_m2_s']:g} | {tau:.3f} | {r['status']} | {r['steps']} | {float(r['physical_time_s'])*1000:g} | {r['max_speed_over_run']:.6g} | {float(d['min_phi']):.6g} ~ {float(d['max_phi']):.6g} |")
lines += ['', '## 비교 해석','']
for angle in (30,90):
    a=cases[f'theta{angle}_gas_baseline'];b=cases[f'theta{angle}_gas_matched']
    if not a['full_10000_pass'] and b['full_10000_pass']:
        lines.append(f'- {angle}°에서 기존 기체 점성은 실패했지만 기체 점성을 높인 대조군은 10000스텝 통과했습니다. 기체 점성/완화의 민감성을 지지합니다.')
    elif not a['full_10000_pass'] and not b['full_10000_pass']:
        lines.append(f'- {angle}°에서는 두 기체 점성 조건 모두 10000스텝을 통과하지 못해 기체 점성 변경만으로 안정화되지 않았습니다.')
    else:lines.append(f'- {angle}°에서 기존/증가 기체 점성 조건의 결과는 표와 같습니다. 통과 구간 밖의 안정성을 보장하지 않습니다.')
for gas in ('baseline','matched'):
    a=cases[f'theta30_gas_{gas}'];b=cases[f'theta90_gas_{gas}']
    lines.append(f"- 기체 {gas} 조건에서 30°는 {a['status']}({a['steps']}스텝), 90°는 {b['status']}({b['steps']}스텝)입니다. 동일 물리 시간과 전체 진단을 함께 비교해야 합니다.")
lines += ['', '기존 기체 점성의 두 접촉각에서 최종 공간장 이상은 모두 공급원 오른쪽의 내부 기공, 약 (2.68~2.71, 0.205~0.235) mm에 집중됐습니다.',
          '각도 변경에 따라 입자 접촉 계면도 달라질 수 있으나, 이번 시험은 각도 변경만으로 기존 기체 점성의 불안정이 사라지지 않았음을 보여줍니다.']
lines += ['', '## 해석 한계와 후속 조치','',
          '- 90°는 중립 젖음 대조군이며 PhaseFieldCurvedWall 자체를 제거하지 않았습니다. 각도에 따른 계면 경로도 변하므로 경계 알고리즘 오류를 단독 확정하지 않습니다.',
          '- 기체 동점성 1e-6 m²/s 대조군은 기체 NS tau도 0.503에서 0.53으로 바꿉니다. 물리적 감쇠와 수치적 완화를 동시에 바꾸므로 충돌 알고리즘만의 효과로 해석할 수 없습니다.',
          '- 기체 밀도는 100 kg/m³로 유지했습니다. 증가 기체 점성 대조군은 실제 물–공기 물성이 아니며 생산 기본값으로 채택하지 않았습니다.',
          '- 기본 접촉각 30°는 모든 고체 입자와 하단 고체 벽에 적용됩니다. 고체의 밀도/탄성/열 물성은 여전히 설정하지 않았고 고체 역학은 풀지 않습니다.',
          '- phi clipping, 속도 중단 기준 변경 또는 phi 허용 범위 확대는 적용하지 않았습니다.',
          '- PASS는 2.5 ms까지의 유한성·phi/밀도/tau 범위 및 저장 압력/VTK 검사입니다. 질량 수지와 물리 정확성, 장시간 안정성까지 검증한 결과가 아닙니다.',
          '- 접촉각은 앞으로 30°를 유지하고, 기체 점성은 원인 분리용 값과 실제 적용할 물성을 구분해 후속 검증해야 합니다.',
          '- 같은 30°에서 고정된 실제 기체 물성을 유지한 충돌 모델/소스 처리 대조, 정적 계면의 기생 유동·접촉각·Laplace 압력 및 상장 질량 수지 검증을 우선합니다.',
          '', '각 실행 명령/소요 시간: experiment_summary.json. 비교: comparison.csv. 소스 스냅샷: source_at_execution/. 재생성: python3 analyze_wetting_gas.py.']
(BASE/'REPORT.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(checks,indent=2))

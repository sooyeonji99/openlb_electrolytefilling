"""Audit and report the completed real-air 30000-step validation."""
import csv,hashlib,json,math,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent
BASE=ROOT/'results/air_30k_v1'
SELECTED=BASE/'air_corrected_rlb'


def csvrows(path):
    with path.open() as f:return list(csv.DictReader(f))

if __name__=='__main__':
    result=csvrows(SELECTED/'multiphase_validation.csv')[0]
    cfg=dict(csv.reader((SELECTED/'run_config.csv').open()))
    diagnostics=csvrows(SELECTED/'multiphase_diagnostics.csv')
    assert result['status']=='PASS' and int(result['steps'])==30000
    assert len(diagnostics)==30001 and all(int(r['step'])==i and r['finite']=='1' for i,r in enumerate(diagnostics))
    required={'rho_gas_kg_m3':1.2,'nu_gas_m2_s':1.5e-5,'rho_liquid_kg_m3':1000,'nu_liquid_m2_s':1e-6,'contact_angle_deg':30,'pressure_pa':2500,'pressure_ramp_s':.001,'phase_tau':.53,'tau_liquid':.53,'tau_gas':.95,'dx_m':5e-6,'dt_s':2.5e-7,'sigma_N_m':1e-5}
    for key,value in required.items():assert math.isclose(float(cfg[key]),value,rel_tol=1e-12,abs_tol=1e-18),(key,cfg[key])
    assert cfg['ns_collision']=='crlb' and cfg['phase_projection']=='conservative' and cfg['material_interpolation']=='bounded'
    assert float(result['reservoir_pressure_pa'])==2500
    assert (SELECTED/'circles.csv').read_bytes()==(ROOT/'results/geometry_dx05/circles.csv').read_bytes()
    for name in ['geometry_validation.csv','particle_validation.csv']:assert csvrows(SELECTED/name)[0]['status']=='PASS'
    # Keep field verification independent of solver's own PASS flag.
    with (BASE/'air_corrected_rlb_vtk_check.log').open('w') as f:
        process=subprocess.run(['python3',str(ROOT/'verify_vtk.py'),str(SELECTED)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    assert process.returncode==0,'VTK pressure/finite check failed'
    vtks=csvrows(SELECTED/'vtk_validation.csv');assert len(vtks)==301 and all(r['finite_fluid_fields']=='1' for r in vtks)
    summary={'case':SELECTED.name,'status':'PASS','steps':30000,'physical_time_s':float(result['physical_time_s']),'target_pressure_pa':2500,'final_pressure_pa':float(result['reservoir_pressure_pa']),'gas_density_kg_m3':1.2,'gas_nu_m2_s':1.5e-5,'gas_mu_Pa_s':1.8e-5,'liquid_tau':.53,'phase_tau':.53,'gas_tau':.95,'dt_s':float(cfg['dt_s']),'phase_mobility_m2_s':(.53-.5)*float(cfg['dx_m'])**2/(3*float(cfg['dt_s'])),'max_speed_over_run':max(float(r['max_speed_lattice']) for r in diagnostics),'raw_min_phi':min(float(r['raw_min_phi']) for r in diagnostics),'raw_max_phi':max(float(r['raw_max_phi']) for r in diagnostics),'final_min_phi':float(diagnostics[-1]['min_phi']),'final_max_phi':float(diagnostics[-1]['max_phi']),'projection_l1_total':sum(float(r['phase_projection_l1']) for r in diagnostics),'projection_l1_max_step':max(float(r['phase_projection_l1']) for r in diagnostics),'projection_mass_residual_max':max(abs(float(r['phase_projection_mass_residual'])) for r in diagnostics),'diagnostic_rows':len(diagnostics),'vtk_snapshots_checked':len(vtks),'max_speed_saved_vtk':max(float(v['max_fluid_speed_lattice_vtk']) for v in vtks),'geometry_identical':True,'no_temperature_evolution':True}
    assert summary['projection_mass_residual_max']<1e-6
    summary['low_mach_speed_screen_pass']=summary['max_speed_over_run']<=.1
    (BASE/'accepted_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    comparisons=[]
    for folder in sorted(BASE.iterdir()):
        if not folder.is_dir() or not (folder/'run_config.csv').exists() or (folder/'EXCLUDED.txt').exists():continue
        c=dict(csv.reader((folder/'run_config.csv').open()));d=csvrows(folder/'multiphase_diagnostics.csv');last=d[-1]
        status='RUNNING';step=int(last['step'])
        final=folder/'multiphase_validation.csv'
        if final.exists():status=csvrows(final)[0]['status']
        else:
            log=folder.with_suffix('.log').read_text() if folder.with_suffix('.log').exists() else ''
            if 'Phase mass incompatible' in log or 'Nonfinite phase' in log:status='FAIL_EXCEPTION'
        comparisons.append({'case':folder.name,'status':status,'last_logged_step':step,'requested_steps':c['requested_steps'],'collision':c['ns_collision'],'source_mode':c['source_mode'],'projection':c.get('phase_projection','off'),'material_interpolation':c.get('material_interpolation','linear'),'liquid_tau':c['tau_liquid'],'phase_tau':c['phase_tau'],'dt_s':c['dt_s'],'target_pressure_pa':c['pressure_pa'],'pressure_ramp_s':c['pressure_ramp_s'],'max_u_over_run':max(float(r['max_speed_lattice']) for r in d)})
    with (BASE/'comparison.csv').open('w') as f:
        w=csv.DictWriter(f,list(comparisons[0]));w.writeheader();w.writerows(comparisons)
    (BASE/'comparison.json').write_text(json.dumps(comparisons,indent=2)+'\n')
    lines=['# 상온 공기 물성 · 30000스텝 검증','',
      '공기 대표 상수: rho=1.2 kg/m³, nu=1.5e-5 m²/s, mu=1.8e-5 Pa·s. 온도 변화는 풀지 않습니다.',
      '상온·대기압 기준의 비압축성 참조 물성 근사이며, 압축성 공기 EOS/에너지 방정식을 추가한 것은 아닙니다. 압력 0 Pa는 외곽의 기준 게이지 압력입니다.',
      '액체 rho=1000 kg/m³, nu=1e-6 m²/s, 기존 검증용 계면장력 1e-5 N/m, 접촉각 30°, 5 µm 확대 기공 형상을 유지합니다. 계면장력을 실제 물–공기 값으로 교체하지 않았습니다.',
      '', '## 완료 조건','',
      '선택 실행: air_corrected_rlb. 2500 Pa로 1 ms 선형 ramp 후 유지, 보정 RLB + 보존형 phi 재분배 + 물성 환산 범위 보호.',
      'NS 액체 tau=0.53, 상장 tau=0.53, 기체 tau=0.95. 원래 dt=0.25 µs와 상장 mobility=1e-6 m²/s를 유지했습니다.',
      f"**30000스텝 = {summary['physical_time_s']*1000:g} ms 완료**, 최종 공급 압력 2500 Pa. 전 스텝 진단 {summary['diagnostic_rows']}행 및 VTK {summary['vtk_snapshots_checked']}개 검사 통과.",
      f"전체 구간 NS 최대 격자 속도: {summary['max_speed_over_run']:.8g}, 저장 VTK 최대 격자 속도: {summary['max_speed_saved_vtk']:.8g}. 최종 phi: [{summary['final_min_phi']:.8g}, {summary['final_max_phi']:.8g}].",
      '', '## 필요한 수치 변경','',
      '- 기존 선형 밀도 혼합에서는 공기 밀도비에서 작은 phi 음수 오차만으로 밀도가 음수가 되어 4스텝에 중단됐습니다.',
      '- 밀도/점성 환산용 phi와 그 밀도 기울기를 물리 범위에 맞추는 보호를 추가했습니다. 이 보호만으로는 안정화되지 않았습니다.',
      '- 상장 분포의 범위 보정은 source를 제외한 pore+출구의 보정 전 총 phi를 보존하면서 diffuse-interface 가중치로 재분배합니다. 분포함수에 w_i*delta_phi를 더해 1차 모멘트를 유지합니다. 원격 계면 사이 재분배가 가능하므로 물리적 정확성은 별도 검증 대상입니다.',
      '- OpenLB risingBubble2d가 사용하는 CorrectedMultiPhaseIncompressibleRLBdynamics/CorrectedLiangPostProcessor를 바탕으로 보정 모델을 적용했습니다. 밀도 Laplacian과 보정된 힘의 밀도 분모, native phase correction을 일관되게 사용합니다. 국소 파일 physical_range_coupling.h에 구현했고 원본 OpenLB 라이브러리는 수정하지 않았습니다.',
      '- 단순 clipping과 달리 총 phi를 보존합니다. 원시 phi 범위, 보정 L1, 보정 총량 잔차를 매 스텝 기록합니다. 총량이 보존 가능한 범위를 벗어나면 예외로 중단합니다.',
      f"보정 전 raw phi 범위: [{summary['raw_min_phi']:.8g}, {summary['raw_max_phi']:.8g}]. 최대 스텝 보정 L1: {summary['projection_l1_max_step']:.8g}, 누적 L1: {summary['projection_l1_total']:.8g} (phi×cell 단위).",
      f"보정 총량 잔차 최대: {summary['projection_mass_residual_max']:.8g}. 이 값은 재분배 연산의 보존 오차이며 공급/배출을 포함한 전체 질량 수지 검증과는 다릅니다.",
      '', '## 제외/실패 후보','',
      'comparison.csv/json에 모든 유효 후보를 기록했습니다. 낮은 dt로 mobility가 바뀐 시험은 최종 조건으로 채택하지 않았습니다.',
      '일부 보정 시험은 보정 후 phi가 범위 내여도 속도와 원시 phi가 발산해 질량 제약 예외로 종료됐습니다. 마지막 CSV의 finite=1만으로 PASS로 분류하지 않았습니다.',
      'prebuild_aborted_old_parameters는 빌드 완료 전에 잘못 시작된 이전 물성 시험으로 중단·제외했으며 공기 검증 결과가 아닙니다.',
      '', '## 검증 범위','',
      '이번 PASS는 7.5 ms 동안의 수치 진행 및 입력/형상/저장 필드 검사입니다. 실제 물리 정확성, 접촉각 재현, Laplace 압력, 상장 mobility/계면 두께/격자 수렴과 전체 질량 수지는 후속 검증해야 합니다.',
      '고체는 고정 no-slip·불투과 벽이며 고체 밀도/탄성/열 물성은 여전히 미설정입니다.',
      '실제 실행 당시 소스: source_corrected/. 실행 명령은 재현 스크립트 및 run_config.csv를 확인합니다. 분석 재생성: python3 analyze_air_validation.py.']
    (BASE/'REPORT.md').write_text('\n'.join(lines)+'\n')
    (BASE/'completion.json').write_text(json.dumps({'complete':True,'selected_case':SELECTED.name,'steps':30000,'physical_time_s':summary['physical_time_s'],'vtk_snapshots_checked':len(vtks),'status':'PASS'},indent=2)+'\n')
    print(json.dumps(summary,indent=2))

# Cell_horizontal

4.3 mm(가로 x) × 5 mm(세로 y)의 2D 직사각형 electrolyte filling 모델입니다.
하단 중앙 x=1.65–2.65 mm, y=0–1 mm를 지속 공급 reservoir로 지정합니다.
수학적 point source는 격자와 공급 유량에 의존하는 특이점이므로 1 mm × 1 mm
영역을 사용했습니다. 외곽 bounce-back 벽 두께 때문에 하단 약 2dx는 벽이고,
reservoir는 그 바로 위부터 시작합니다. 공급은 내부 압력/phi bath로 구현하며
매 스텝 pressure 및 phi=1을 재설정합니다. NS 공급원은 정지 reservoir의
지정 압력 평형 분포로 갱신하고 CH는 phi=1로 갱신합니다. 고정 유량 공급은 아닙니다.

## 유지한 기준 조건

직전 mm 규모 연구 `NMC_timecost/multiphase_smoke.cpp`의 conservative phase-field,
Liang coupling, incompressible BGK와 닫힌 외곽 벽을 기준으로 했습니다.
`NMC_physic_added`의 LP30 프로필/중력/10 Hz 진동은 이 기준 모델에 포함되지
않으므로 여기에도 적용하지 않았습니다. 이전 연구가 LP30 전체 구동 조건을
의미한다면 그 조건과 동일한 비교라고 해석하면 안 됩니다.

| 조건 | 값 |
|---|---:|
| 액체/기체 밀도 | 1000 / 100 kg/m³ |
| 액체/기체 동점도 | 1e-6 / 1e-7 m²/s |
| 액체/기체 tau | 0.53 / 0.503 (기존은 0.8 / 0.53) |
| 표면장력 | 1e-5 N/m |
| 액체 목표 접촉각 | 60° |
| 계면 폭 | 4dx |
| 목표 다공성 영역 공극률 | 0.30 |
| 고체 배치 pitch / 반경 | 50 µm / 약 23.60 µm |
| 난수 seed | 20261014 |
| 기본 격자 | 5 µm |
| 지속 공급 gauge 압력 | 1000 Pa (시험값) |

원형 고체 위치는 각 50 µm 타일에서 x/y 독립 난수로 ±0.75 µm 이동합니다.
겹침과 기공 연결 손실을 피하기 위한 **제한된 무작위 배열**이며 완전히 무작위인
입자 packing이나 기존 submicron NMC 입도분포의 재현은 아닙니다. 모든 고체는
동일 반경입니다. 공급원을 침범하는 입자는 원 전체를 제거하며 원을 잘라내지
않습니다. `circles.csv`에 좌표와 제거 여부를 기록합니다.

공급원은 Pa를 converter의 압력 변환계수로 나눈 값을 `ns.iniEquilibrium`에
전달합니다. 해당 다상 모델에서 일반 `defineRho` 기반 setter는 목표 압력과
실제 압력이 다를 수 있어 사용하지 않습니다. `run_config.csv`에는 물리 압력과 격자 압력/변환계수를
기록하고 VTK에는 실제 계산 pressure_phys_Pa를 출력합니다.

## 빌드와 실행

설치된 OpenLB의 config.mk 설정(MPI 포함)을 상속하며 루트 설정은 수정하지 않습니다.

```bash
make -j2
./geometry_validate 2 results/geometry_dx02
./geometry_validate 5 results/geometry_dx05
./main --tau 0.53 --reservoir-pressure-pa 1000 --steps 200 \
  --diagnostic-every 1 --output-every 100 --output results/verified_1000Pa
./main --tau 0.53 --reservoir-pressure-pa 0 --steps 200 \
  --diagnostic-every 1 --output-every 100 --output results/verified_0Pa
python3 verify_vtk.py results/verified_1000Pa
python3 verify_vtk.py results/verified_0Pa
python3 analyze.py
```

`--tau`로 시간 간격을 바꿀 때 기체 tau도 함께 재계산해 두 상의 물리
동점도를 유지합니다. tau가 0.5에 너무 가까우면 별도 안정성 검증이 필요합니다.

MPI 환경에서는 `mpirun -np 2 ./main ...`으로 실행할 수 있습니다.
모든 local block에 동일한 물리 좌표 형상 함수를 적용합니다.
출력 폴더를 매 실행 다르게 지정해 결과를 보관하세요.

## 출력과 검증 범위

- `geometry.vtk`: 0 외부, 1 pore, 2 solid/외곽 벽, 3 공급 reservoir.
- `geometry_validation.csv`: 공극률 및 reservoir와 연결된 기공/상단 도달 검사.
- `circles.csv`, `run_config.csv`: 난수 입자 목록 및 계산 설정.
- `multiphase_diagnostics.csv`: 시간, 기공 평균 phi, phi 범위, 최대 격자 속도.
- `multiphase_validation.csv`: 유한성/phi [-0.1,1.1]/속도 <0.1 검사.
- `vtkData/cell_horizontal.pvd`: ParaView 시간별 material/phi/유속/압력.

ParaView에서 material=1 또는 3을 선택한 뒤 phi_fluid_only를 0–1로 표시하세요.
고체/외부 phi는 의도적으로 NaN으로 마스킹합니다. 기공 평균 phi는 material 1만
포함하며 reservoir의 초기 액체량은 충전율에 포함하지 않습니다. VTK의 시간은
스텝이며 물리 시간은 CSV에서 확인하세요.

닫힌 외곽에는 gas vent가 없습니다. incompressible 기체 모델은 실제 기체의
열역학적 압축을 재현하지 않습니다. reservoir의 강제 재설정으로 액체가 추가될
수 있으므로 충전율 변화만으로 정확한 보존 수지를 검증했다고 해석하지 않습니다.
PASS는 짧은 통합 안정성 검사이며 실제 LP30/NMC 침투나 장시간 충전 검증이 아닙니다.

기존 tau=0.8의 5 µm 격자는 약 86만 노드, dt=2.5 µs입니다.
고압 계산에 사용하는 기본 tau=0.53에서는 dt=0.25 µs이고 10초에 약
4000만 스텝이 필요합니다. 2 µm에서는 같은 tau로 dt=0.04 µs, 약 2억5000만
스텝이 필요합니다. 전체 장시간 계산은 이번 실행에 포함하지 않았습니다.
10 µm 형상은 공극률 오차가 허용범위(±0.04)를 넘어서 FAIL이므로 계산에 사용하지 않습니다.

## 실행 결과 (2026-10-14)

사용자가 지정한 비교 기준은 NMC_timecost의 검증용 물성·닫힌 외곽입니다.
물성과 형상 크기는 그대로 두고 고압 안정성을 위해 dt만 축소했습니다.

| 격자 | 공극률 | 공급원 연결 기공 비율 | 형상 검사 |
|---|---:|---:|---|
| 2 µm | 0.296452 | 0.999958 | PASS |
| 5 µm | 0.307442 | 1.000000 | PASS |
| 10 µm | 0.358459 | 1.000000 | FAIL (공극률 오차) |

5 µm / dt=0.25 µs / 200스텝 = 50 µs의 가압/무가압 비교입니다.
두 상의 동점도는 기존 모델과 동일합니다.

| 지속 공급 압력 | 초기 기공 액상 분율 | 최종 기공 액상 분율 | 최대 격자 속도 | 검사 |
|---|---:|---:|---:|---|
| 0 Pa | 0.0004293543 | 0.0014119543 | 0.000025262 | PASS |
| 1000 Pa | 0.0004293543 | 0.0014223265 | 0.007964867 | PASS |

두 조건 모두 매 스텝 진단을 저장했습니다. 각 3개 VTK의 유체 phi/압력/유속
유한성을 확인했습니다. 실제 가압 공급원 압력은 약 999.99976–1000 Pa입니다.
가압 조건 최종 액상 분율은 무가압보다 약 0.73% 높지만 차이가 작고,
두 조건에 공통된 초기 계면 완화가 포함되어 있습니다. 이 짧은 결과로
장시간 침투율 또는 전체 충전 시간을 확정할 수 없습니다.

최종 결과는 `results/verified_1000Pa/`, `results/verified_0Pa/`에 있습니다.
ParaView는 각 폴더의 `vtkData/cell_horizontal.pvd`를 열면 됩니다.
`geometry_preview.svg`는 이상적인 외곽/입자/source 배치 미리보기입니다.
정확한 래스터 형상은 geometry.vtk와 OpenLB material 출력으로 확인하세요.

`results/superseded*`는 압력 단위/공급원 setter를 점검하던 중간 실행입니다.
실제 1000 Pa 검증 결과로 사용하지 마세요. 최종 분석 스크립트는 verified
실행만 집계합니다. 고압·기존 dt 조합의 속도 제한 초과는 이 중간 점검에서
확인했으며 해당 setter 출력과 지정 압력이 일치하지 않아 정량 비교에서 제외했습니다.

최종 solver의 기본 옵션 실행(20스텝)도 PASS이며 1000 Pa와 dt=0.25 µs를
확인했습니다. 실제 OpenLB material과 생성 격자의 861,861개 셀을 대조한
불일치는 0입니다. 최종 구현으로 tau=0.8 / 1000 Pa를 다시 시험한
`results/stability_tau08/`는 2스텝에서 최대 격자 속도 0.291852로 FAIL입니다.
따라서 기본 tau=0.53을 사용합니다. 이는 물성이 바뀐 것이 아니라 시간
간격과 두 상의 완화시간을 함께 줄인 것입니다.

## 닫힌 셀 침투 정지 / 잔류 기체 검증

1000 Pa 공급·닫힌 외곽·기존 물성·dx=5 µm·dt=0.25 µs를 유지하고
0.5 ms(2000스텝)까지 연장하는 검증을 추가했습니다.

```bash
./main --steps 2000 --diagnostic-every 1 --output-every 500 \
  --output results/closed_cell_0p5ms_final > results/closed_cell_0p5ms_final.log 2>&1
python3 analyze_closed_cell.py results/closed_cell_0p5ms_final
python3 analyze.py
```

`analyze_closed_cell.py`는 다음 파일을 생성합니다.

- `closed_cell_history.csv`: 기공 평균 phi 및 raw 기체 면적/비율 추정.
- `closed_cell_analysis.json`: 안정성, 침투 정지 후보 판정, 마지막 VTK 시각.
- `residual_gas_components.csv`: 마지막 직렬 VTK에서 phi<0.5 기체 영역의
  4방향 연결 성분, 면적, 중심 좌표, bounding box.

침투 정지 후보는 마지막 0.1 ms 구간 3개에서 각 충전율 범위가 **절대 1e-6
및 해당 충전율의 0.1% 이하**를 모두 만족할 때만 표시합니다. 각 구간에
3개 이상 진단 표본과 시간 길이의 95% 이상이 필요합니다. 이 기준은 임시
판정 기준이며 사용자가 `--window-s`, `--absolute-tolerance`,
`--relative-tolerance`로 지정할 수 있습니다.

`NUMERICAL_FAILURE`는 수치 검사 실패, `NOT_YET_DEMONSTRATED`는 시간 부족
또는 변화가 계속됨, `PLATEAU_CANDIDATE`는 추가 확인이 필요한 정지 후보입니다.
실패 상태에서는 plateau로 판정하지 않습니다. 최종 VTK가 종료 시각과 다르면
`is_final_step=false`를 표시합니다. MPI 분할 VTK의 기체 연결 성분은 halo를
합치지 않고 계산하면 잘못되므로 이 분석은 직렬 VTK만 지원합니다.

중단 시점에도 VTK를 저장하도록 solver를 보완했습니다. 잔류 기체 면적은
2D 상장 proxy이며 기체 질량이 아닙니다. phi<0.5 성분 수만으로 물리적인
기포 포획/압축을 검증했다고 해석하지 않습니다. 닫힌 셀 모델에도 내부
reservoir 강제 갱신을 통한 공급은 존재하며, 기체는 비압축성입니다.

초기 연장 실행 `results/closed_cell_0p5ms/`는 430스텝(0.1075 ms)에서
최대 격자 속도 0.14384로 중단했습니다. 충전율은 계속 증가 중이었으며,
**닫힌 셀 침투 정지는 검증되지 않았습니다**. 최종 snapshot 저장 기능을
추가한 재실행 결과는 아래에 기록합니다. 이 실패를 해결하기 전에 5/50 ms
계산으로 단순 연장하지 않습니다.

### 중단 시점까지 재검증한 결과

`results/closed_cell_0p5ms_final/`에서 매 스텝 검사한 결과,
**415스텝 = 0.10375 ms**에 최대 격자 속도 **0.108955 > 0.1**로 중단했습니다.
목표 0.5 ms에는 도달하지 못했습니다. 초기 실행은 10스텝마다 검사했으므로
430스텝에서 실패를 발견했고, 이번에는 더 이른 최초 기준 초과를 확인했습니다.

- 기공 액상 분율: 0.000429354 → 0.001704376 (약 0.17044%).
- raw 기체 면적 비율: 약 99.82956% (기공 면적 기준).
- phi<0.5 기체 노드: 249,906개, 4방향 연결 성분 1개.
- 공급원과 연결된 phi>=0.5 액체 기공 노드: 484개.
- 종료 시각 VTK: `vtkData/cell_horizontal.pvd`, step 415 포함.
- 초기/종료 VTK의 유체 필드는 모두 유한하며 공급원 압력은 1000 Pa입니다.

**비유한 필드가 발생한 것이 아니라 격자 속도 안전 기준을 초과한 것입니다.**
침투 정지나 물리적 잔류 기포 상태에 도달했다는 결과가 아닙니다. 최종 기체
연결 성분은 대부분 아직 채워지지 않은 기공 영역이며 기포 포획 증거가 아닙니다.
분석 상태는 `NUMERICAL_FAILURE`입니다. 장시간 정지 검증 전에 고압 구동에서
시간 간격/충돌 모델 및 공급원 처리를 검토하고 수치 기준을 만족해야 합니다.

분석 판정은 임시 fixture로 일정 충전율/계속 증가/수치 실패의 세 경우를
확인했습니다. solver 빌드, Python 문법 검사 및 실제 종료 snapshot 분석을
완료했습니다. 기존 형상·물성·공급 압력·외곽 조건은 유지했습니다.

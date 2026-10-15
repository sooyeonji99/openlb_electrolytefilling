# Cell_horizontal — random mixed particles

현재 코드는 기존 단일 직경 약 47.2 µm 타일 배열을 덮어쓰고, 사용자 지정에
따라 **직경 50–200 µm** 원형 고체 혼합 배열로 변경했습니다. 최초 요청의
5–20 µm 대신 기존보다 큰 50–200 µm를 사용하기로 확인한 조건입니다.

직사각형은 x=0–4.3 mm, y=0–5 mm이며 하단 중앙 x=1.65–2.65 mm,
y=0–1 mm는 electrolyte reservoir입니다. 원형 고체는 reservoir를 침범하지
않으며 외곽은 닫힌 벽입니다. 목표 공극률 0.30은 reservoir와 외곽 벽을
제외한 다공성 영역에 적용합니다.

## 새 입자 배열

- 직경을 50–200 µm 구간에서 균등 난수로 샘플링합니다.
- 초기 중심 좌표를 무작위로 생성하고 반복 이동으로 겹침을 제거합니다.
- 원끼리의 물리 표면 간격은 최소 약 10 µm로 유지합니다.
- 공급원과 원 표면 사이에도 최소 약 10 µm 간격을 유지합니다.
- seed=20261014, 1,035개 입자, 직경 범위 약 50.074–199.983 µm.
- 동일한 물리 입자 목록을 2/5 µm 격자에 래스터화합니다. 격자에 따라
  원 크기나 위치를 다시 생성하지 않습니다.

목표 고체 면적은 닫힌 외곽 부근 10 µm 물리 packing margin과 reservoir를
제외한 영역으로 계산합니다. 최종 래스터 공극률은 노드 샘플링 및 외곽 벽
두께 때문에 정확히 0.30과 같지는 않으며 허용 오차는 절대 ±0.005입니다.
원 사이 간격 10 µm는 물리 원 표면 거리이며 래스터 기공 전체의 최소 throat
폭이나 유동 정확성을 보장하는 값은 아닙니다. 입자-외곽 벽 간격은 별도입니다.

| 격자 | 입자 직경 해상도 | 공극률 | 공급원 연결 기공 비율 | 검사 |
|---|---|---:|---:|---|
| 5 µm | 약 10–40셀 | 0.301618 | 1.000000 | PASS |
| 2 µm | 약 25–100셀 | 0.304211 | 1.000000 | PASS |

모든 입자 쌍을 검사한 최소 표면 간격은 9.999917 µm이며 겹침은 없습니다.
두 격자의 `circles.csv`는 동일합니다. 5 µm보다 큰 격자는 현재 허용하지 않습니다.

## 유지한 물리/계산 조건

`NMC_timecost`의 conservative phase-field / Liang coupling / incompressible
BGK 모델 및 닫힌 외곽을 유지합니다. LP30, 중력, 진동, gas vent는 적용하지
않습니다. 액체/기체 밀도 1000/100 kg/m³, 동점도 1e-6/1e-7 m²/s,
표면장력 1e-5 N/m, 목표 접촉각 60°, 계면 두께 4dx입니다.
기본 dx=5 µm, tau_liquid=0.53, tau_gas=0.503, dt=0.25 µs입니다.
1000 Pa reservoir 압력 평형 분포와 phi=1을 매 스텝 갱신합니다.

입자 크기와 간격을 바꾸어도 이전 속도 제한 초과의 원인이 병목이었다고
확정할 수는 없습니다. 기존 모델은 415스텝에서 속도 기준을 초과했으며,
새 형상의 안정성과 침투 거동은 새 결과로 따로 판단해야 합니다.

## 실행

```bash
make -j2
./geometry_validate 5 results/mixed_geometry_dx05
./geometry_validate 2 results/mixed_geometry_dx02
python3 preview_geometry.py

# 새 형상 단기 가압 시험: 0.125 ms
./main --steps 500 --diagnostic-every 1 --output-every 250 \
  --output results/mixed_pressure_1000Pa_0p125ms \
  > results/mixed_pressure_1000Pa_0p125ms.log 2>&1
python3 verify_vtk.py results/mixed_pressure_1000Pa_0p125ms
python3 analyze_closed_cell.py results/mixed_pressure_1000Pa_0p125ms
python3 analyze.py
```

기본 결과 디렉터리도 `results/mixed_pressure_1000Pa`로 변경했습니다.
과거 `verified_*`, `closed_cell_*`, `superseded*` 폴더는 이전 형상의 결과입니다.
기존 README는 `results/geometry_v1_README.md`에 보관했습니다. 기존 결과를
새 형상의 검증 결과로 해석하지 마세요.

0.2초 목표는 같은 dt에서 `--steps 800000`입니다. 속도 <0.30, phi 범위
[-0.1,1.1], 필드 유한성 검사 실패 시 조기 중단하며 중단 시점 VTK도 저장합니다.
checkpoint 재시작 기능은 없으므로 모든 실행은 초기 상태부터 시작합니다.

## 출력

- `geometry.vtk`: material 1=pore, 2=solid/닫힌 벽, 3=reservoir, 0=외부.
- `circles.csv`: 물리 입자 좌표 및 반경(m); 직경은 반경의 두 배.
- `particle_validation.csv`: 입자 수, 직경 범위, 모든 쌍의 최소 간격 검사.
- `diameter_histogram.csv`: 25 µm 구간별 입자 수.
- `geometry_validation.csv`: 공극률, 연결성 및 상단 내부 도달 검사.
- `run_config.csv`: `geometry_version=random_mixed_disks_v2`와 물성/입자 설정.
- `vtkData/cell_horizontal.pvd`: 시간별 material, phi, 실제 압력, 물리 속도.
- `closed_cell_analysis.json`: 수치 검사 및 침투 정지 후보 판정.

`geometry_preview.svg`는 이상적인 물리 원 표면 미리보기입니다. 실제 solver
래스터 형상은 VTK material로 확인하세요. phi를 표시할 때 material 1/3만
선택하세요. 고체 phi의 NaN은 의도적인 마스킹입니다.

침투 정지 후보는 0.1 ms 구간 3개에서 충전율 범위가 절대 1e-6와 상대 0.1%
기준을 모두 만족할 때만 판정합니다. 수치 검사 실패는 `NUMERICAL_FAILURE`,
시간 부족 또는 지속 변화는 `NOT_YET_DEMONSTRATED`입니다. 잔류 기체는
phi<0.5의 4방향 연결 성분으로 분석하며 2D 면적 proxy입니다. 비압축성
기체·강제 reservoir 모델에서는 이 결과만으로 실제 기체 압축이나 물리적
포획을 검증할 수 없습니다. MPI 분할 VTK의 연결 성분 분석은 지원하지 않습니다.

## 새 형상 가압 시험 결과 (2026-10-14)

`results/mixed_pressure_1000Pa_0p125ms/`에서 1000 Pa 공급으로
500스텝(0.125 ms)을 실행했습니다. 매 스텝의 수치 검사는 **PASS**입니다.
전체 기간 최대 격자 속도는 **0.034956**, 최종 기공 액상 분율은 **0.005097**
(약 0.5097%)입니다. 초기/중간/종료 VTK 유체 필드는 모두 유한하며,
공급원 압력은 약 999.99994–1000 Pa입니다. 실제 OpenLB 재료와 형상 생성
격자의 861,861개 노드 대조에서 불일치는 0입니다.

이는 새 형상 적용과 초기 가압 거동의 확인입니다. 0.2초 안정성이나 병목이
이전 중단의 유일한 원인이었다는 사실을 검증한 결과는 아닙니다. 침투는
계속 진행 중이며 정지 판정은 `NOT_YET_DEMONSTRATED`입니다.

격자 속도 중단 기준을 0.10에서 0.15로 완화했습니다. 과거 결과는 당시
0.10 기준이며 새 실행은 run_config.csv에 0.15 기준을 기록합니다.
phi 범위와 NaN/Inf 중단 조건은 유지합니다.

사용자 지정에 따라 현재 격자 속도 중단 기준은 0.30입니다(격자 마하수
약 0.52). 새 실행은 run_config.csv에 0.30을 기록하며 기존 0.15 기록은
당시 설정입니다. phi 범위와 NaN/Inf 중단 조건은 유지합니다.

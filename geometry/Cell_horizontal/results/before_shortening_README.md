# 압력 ramp 적용

기본 공급 압력은 0 Pa에서 1 ms 동안 **선형으로 10 kPa**까지 증가한 후
유지합니다. `--pressure-ramp-s`로 상승 시간을 바꾸며 0이면 즉시 가압합니다.
하단 벽, 상단·좌우 0 Pa 출구, 입자 형상 및 물성은 유지합니다.
속도 크기 중단은 비활성화 상태이며 NaN/Inf와 phi 범위 검사는 유지합니다.

```bash
./main --outer-boundary open --reservoir-pressure-pa 10000 \
  --pressure-ramp-s 0.001 \
  --steps 30000 --diagnostic-every 1 --output-every 2000 \
  --output results/open_bottomwall_10kPa_ramp1ms_30000 \
  > results/open_bottomwall_10kPa_ramp1ms_30000.log 2>&1

python3 verify_vtk.py results/open_bottomwall_10kPa_ramp1ms_30000
python3 analyze.py
```

현재 dt=0.25 µs에서 4000스텝에 10 kPa에 도달하고, 30000스텝의 총 물리
시간은 7.5 ms입니다. 더 완만한 상승은 `--pressure-ramp-s 0.005`로
지정하면 됩니다(5 ms, 20000스텝). 초기 상태부터 계산하며 이전 계산의
checkpoint를 이어받지는 않습니다.

CSV의 `reservoir_pressure_pa`는 각 행 시각의 적용 압력이고 run_config의
`pressure_pa`/`pressure_lattice`는 최종 목표값입니다. 각 적분 구간은 시작
시각의 압력으로 진행하고 끝에서 새 시각 압력을 적용합니다. 초기 압력은
ramp>0이면 0 Pa입니다. ramp는 선형이며 시간 간격은 변경하지 않습니다.
`verify_vtk.py`도 각 snapshot 시각의 ramp 압력과 실제 공급원 압력을 비교합니다.
ramp가 10 kPa 장시간 안정성을 보장하지는 않습니다.

---

하단 벽 복원: 10 kPa / 30,000스텝 실행은 다음 명령을 사용합니다.

```bash
./main --outer-boundary open --reservoir-pressure-pa 10000 \
  --steps 30000 --diagnostic-every 1 --output-every 2000 \
  --output results/open_bottomwall_10kPa_30000 \
  > results/open_bottomwall_10kPa_30000.log 2>&1
```

현재 `open`은 상단·좌우 개방과 하단 벽을 뜻합니다. 과거 전체 외곽 개방
실행과 구분하도록 `run_config.csv`에 `bottom_boundary=wall`을 기록합니다.
입자 목록과 내부 reservoir 및 물성은 그대로이며 압력 ramp는 추가하지 않았습니다.
기존 10 kPa·전체 외곽 개방 실행은 104스텝에 NaN이 발생했습니다.

# 현재 기본 조건: 열린 외곽

상·좌·우는 0 Pa gauge 압력 출구이며 하단 전체는 벽입니다. 하단 중앙 내부의
1 mm × 1 mm / 1000 Pa electrolyte 공급원은 유지합니다. 입자 좌표·반경,
내부 pore/solid 분포, 물성, dx, dt는 유지합니다. 현재 속도 크기에 따른
중단은 비활성화되어 있으며 속도 NaN/Inf와 phi 범위 검사는 유지합니다.
아래의 과거 닫힌 조건 결과는 열린 경계의 장시간 검증 결과가 아닙니다.

```bash
make -j2
mkdir -p results
./main --outer-boundary open --tau 0.53 --reservoir-pressure-pa 1000 \
  --steps 50000 --diagnostic-every 1 --output-every 2000 \
  --output results/open_all_50000 \
  > results/open_all_50000.log 2>&1

# 종료 후
python3 verify_vtk.py results/open_all_50000
python3 analyze.py
```

50,000스텝은 물리 시간 12.5 ms입니다. 400,000스텝(0.1초)을 목표로
한다면 `--steps 400000`과 새 출력 경로를 사용하세요. 기본 `--outer-boundary`
값은 open이며 닫힌 비교 계산은 `--outer-boundary closed`입니다.
기존 중단 상태를 이어받지 않고 초기 상태부터 계산합니다.

NS에는 검증 예제의 IncompressibleZouHePressure, 상장에는 convective
phase-field 출구를 사용합니다. phi=0을 강제로 고정하는 gas bath가 아니므로
액체도 외곽으로 빠져나갈 수 있습니다. 측면과 상단으로 액체가 빠져나갈 수 있으며 하단을 통한 배출은 차단합니다. 모든 유체/공급원/출구에서 phi 범위와
유한성을 검사하고 기공 평균 충전율은 material 1만 대상으로 계산합니다.

경계는 원래 벽의 안쪽 격자층에 배치합니다. 가장 바깥 한 층은 exterior
(0), 열린 면은 material 4=좌, 5=우, 7=상이며 하단은 material 2입니다.
material 6은 과거 하단 개방 실행에서만 사용했습니다. 네 개의 단일
모서리 격자점과 하단 전체를 material 2로 두어 평면 출구
법선을 정의합니다. 그 외의 pore/solid 및 reservoir 셀은 닫힌 조건과 같습니다.

`verify_vtk.py`는 열린 출구의 압력, 순 대류 유출, 액상 대류 유출 및 역류
노드 수를 snapshot별 `vtk_validation.csv`에 기록합니다. 이 값은 단위 깊이의
2D 유량이며 정확한 상장 보존 수지는 아닙니다. MPI에서는 각 블록 결과이므로
전체 유량으로 합칠 때 겹침 셀 처리가 추가로 필요합니다.
대류 출구는 별도의 외부 유입 상 조성을 지정하지 않습니다. 역류가 지속되면
출구에서의 상 처리에 대한 추가 검증이 필요합니다.

입자 목록은 기존 mixed_geometry_dx05와 바이트 단위로 동일하고 공극률은
0.301618, 공급원 연결 기공 비율은 1.0입니다. 초기 100스텝 열린 계산에서
출구 3,507개, 0 Pa 압력 유지, 순유출을 확인했습니다. 일부 역류도 관측했습니다.
열린 경계가 22,000스텝 부근의 중단을 해결했는지는 아직 검증하지 않았습니다.

`analyze_closed_cell.py`의 닫힌 셀 정지 판정은 열린 계산에 적용하지 않습니다.
열린 조건에서는 배출 액체량도 포함해 충전과 배출을 구분해야 합니다.


## 물리 조건 및 형상

직사각형 4.3 × 5 mm, 하단 중앙 reservoir 1 × 1 mm, 원형 입자 1,035개,
직경 약 50.074–199.983 µm, 최소 입자 간격 약 10 µm, seed=20261014입니다.
5 µm 격자의 내부 공극률은 0.301618, 공급원 연결 기공 비율은 1.0입니다.
2 µm 격자에서도 동일한 물리 입자 목록을 사용하며 공극률은 0.304211입니다.
열린 경계로 바꿀 때 내부 pore/solid/reservoir 변경 셀은 0개입니다.

NMC_timecost 기준 검증용 물성: 액체/기체 밀도 1000/100 kg/m³,
동점도 1e-6/1e-7 m²/s, sigma=1e-5 N/m, 목표 접촉각 60°, 계면 두께 4dx.
LP30, 중력, 진동, 실제 기체 압축/용해 모델은 포함하지 않습니다.
기본 dx=5 µm, dt=0.25 µs, tau_liquid=0.53, tau_gas=0.503입니다.

`run_config.csv`의 outer_boundary로 열린/닫힌 실행을 구분합니다.
`circles.csv`는 좌표와 반경(m), `particle_validation.csv`는 입자 간격 검사,
`geometry_validation.csv`는 내부 공극률 및 연결성 검사입니다.
`geometry_preview.svg`는 이상적인 물리 원 표면을 표시하며 실제 raster는
`geometry.vtk` 및 solver material 출력으로 확인합니다.

## 비교 검증 및 결과 보관

```bash
./geometry_validate 5 results/open_geometry_dx05
./main --outer-boundary open --steps 500 --diagnostic-every 1 \
  --output-every 250 --output results/open_all_500steps \
  > results/open_all_500steps.log 2>&1
python3 verify_vtk.py results/open_all_500steps
```

기존 `mixed_*`, `verified_*`, `closed_cell_*` 결과는 해당 폴더의 당시
설정으로 해석해야 합니다. 일부 이전 경로에는 재실행으로 파일이 혼재할 수
있으므로 새 실행에는 반드시 고유한 결과 경로를 사용하세요.
닫힌 조건 설명은 `results/geometry_v2_closed_README.md`에 보관했습니다.
이전 닫힌 실행은 0.30 기준에서도 22,046스텝(5.5115 ms)에 속도 제한을
초과했습니다. 열린 경계에서 그 시각 이후의 안정성은 별도로 확인해야 합니다.

## 열린 경계 적용 검증 결과 (2026-10-15)

`results/open_all_500steps/`에서 500스텝(0.125 ms)을 실행해 수치 검사를
통과했습니다. 최대 격자 속도는 약 0.03031, 최종 기공 평균 phi는
0.00516738입니다. 3개 snapshot에서 유체 필드는 유한하고 출구 압력은
0 Pa로 유지됐으며 공급원 압력은 약 999.99994–1000 Pa였습니다.
최종 외곽 순 대류 유출은 약 1.689e-5 m²/s이며 역류 노드는 424개였습니다.
이는 초기 경계 동작 검사이며 장시간 역류·상장 보존·안정성 검증은 아닙니다.

`--outer-boundary closed` 20스텝 회귀에서는 이전 동일 형상 계산의 기공
평균 phi와 최대 속도가 모든 비교 스텝에서 정확히 일치했습니다.
입자 목록도 동일하며 boundary_comparison.json에 내부 변경 0개를 기록했습니다.

사용자 지정에 따라 속도 크기에 따른 중단 및 최종 FAIL 조건을 제거했습니다.
최대 속도 기록, NaN/Inf 검사, phi 허용 범위 [-0.1,1.1] 검사는 유지합니다.
새 run_config.csv의 max_speed_lattice_limit은 disabled로 기록됩니다.
PASS는 남아 있는 수치 검사를 통과했다는 의미이며 속도 적합성 검증은 아닙니다.

## 하단 벽 복원 후 10 kPa 시험 (2026-10-15)

`results/open_bottomwall_10kPa_check200/`에서 200스텝을 목표로 시험했으나,
119스텝(29.75 µs)에 phi_min=-0.224701이 허용 하한 -0.1을 벗어나 중단했습니다.
최대 격자 속도 통계는 13.886까지 증가했습니다. 속도 크기 중단은 비활성화된
상태이며 이번 직접 중단 원인은 phi 범위 이탈입니다.

이전 104스텝 NaN은 피했지만 고압 조건의 급격한 속도 증가는 해소되지
않았습니다. 하단 출구만이 불안정의 유일한 원인이라고 볼 수 없습니다.
하단 벽 수정은 유지했으며 압력 ramp, dt, 충돌 모델, 물성은 이번 작업에서
변경하지 않았습니다. 10 kPa 장시간 계산에 앞서 추가 안정화 검증이 필요합니다.

`open_bottomwall_geometry_dx05/boundary_comparison.json`에 하단 출구
656개만 material 6에서 벽 2로 변경됐음을 기록했습니다. 입자 목록은 기존과
동일하고 내부 공극률 0.301618 및 공급원 연결 기공 비율 1.0을 유지했습니다.

Ramp 구현 검증: 1 ms ramp의 초기 200스텝(50 µs, 500 Pa까지)은 PASS입니다.
별도 10 µs 단기 ramp 시험에서는 VTK 공급원 압력이 0→5000→10000 Pa로
전환되고 이후 유지되는 것을 확인했습니다. 모든 기록 시각의 CSV 압력은
해석식과 일치했습니다. 기본 1 ms ramp로 10 kPa에 도달한 뒤 장시간 유지하는
계산은 아직 수행하지 않았습니다.

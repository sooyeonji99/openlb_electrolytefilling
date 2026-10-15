# phi 안정성 검증 결과 추가

현재 실행 조건은 상장 tau를 이전 값 **0.53**으로 사용합니다.
아래 명령은 2 µm, 20,000스텝, 100스텝 간격 출력이며,
tau=1.0 결과와 구분되는 폴더에 저장합니다.

```bash
./main --dx-um 2 --phase-tau 0.53 \
  --outer-boundary open --reservoir-pressure-pa 10000 \
  --pressure-ramp-s 0.001 --steps 20000 \
  --diagnostic-every 1 --output-every 100 \
  --output results/dx02_phase_tau053_20000 \
  > results/dx02_phase_tau053_20000.log 2>&1
```

동일한 축소 형상에서 13개 대조 실행을 완료했습니다.
결과: [안정성 검증 보고서](results/stability_v1/REPORT.md),
[비교 CSV](results/stability_v1/comparison.csv), [검증 옵션](validation/README.md).

기존 기본값은 재현을 위해 유지했습니다. 새 `--phase-tau`, `--ns-collision`,
`--source-mode`, `--sigma`는 검증용 선택 옵션입니다. 상장 tau=1 후보는
5 µm에서 10 kPa 도달 후까지 유한성/phi 검사를 통과했지만 mobility를
바꾸는 조건이므로 물리적으로 검증된 최종값이 아닙니다. 2 µm에서는 같은
후보가 1300스텝까지 phi 검사를 통과했어도 최대 격자 속도가 0.3을 넘어
장시간 안정화로 판정하지 않았습니다. 진단에 min_rho_lattice/min_tau_eff를
추가했고 음의 혼합 밀도와 NS tau<=0.5도 중단 조건으로 기록합니다.
속도 크기 중단 비활성화와 bulk phi 비클리핑은 유지합니다.

---

# Cell_horizontal — 높이 축소

기존 공급원 위 4 mm 구간을 1/3(1.333333 mm)로, 하단 electrolyte 공급원
높이를 1 mm에서 0.5 mm로 축소했습니다. 가로 4.3 mm와 공급원 폭 1 mm는
유지합니다. 명목 총 높이는 1.833333 mm입니다. 공급원 옆 하단 영역에도
기존과 같이 pore/solid가 존재합니다.

격자에 맞춰 전체 높이를 반올림하며 요청 높이와 실제 높이를 run_config의
requested_height_m / height_m에 각각 기록합니다. converter는 가로 폭을
기준으로 설정해 실제 격자 간격과 요청 dx가 일치하도록 했습니다.

| dx | 실제 전체 높이 | 노드 수 (x × y, 경계 포함) | 전체 노드 | 공극률 |
|---|---|---|---:|---:|
| 5 µm | 1.835 mm | 861 × 368 | 316,848 | 0.297851 |
| 2 µm (기본) | 1.834 mm | 2151 × 918 | 1,974,618 | 0.302238 |

두 해상도 모두 공급원 연결 기공 비율 100%, 형상·입자 검사는 PASS입니다.
입자 직경 50–200 µm와 최소 입자 간격 약 10 µm 및 목표 공극률 0.30은
유지하고 새 영역에 맞춰 다시 배치했습니다. 2/5 µm에서 입자 목록은 동일합니다.
이전 긴 영역의 입자 좌표를 수직 압축한 형상은 아닙니다.

기본 공간 해상도를 2 µm로 변경했습니다. 입자 직경은 25–100셀,
최소 입자 간격 약 10 µm는 약 5셀로 표현합니다. 5 µm 비교 계산은
`--dx-um 5`로 지정할 수 있습니다. 기존 861×1001=861,861노드에 비해 노드 수는 약 36.8%입니다.

하단 벽 / 상단·좌우 0 Pa 개방, 검증용 물성, 속도 중단 비활성화,
NaN/Inf 및 phi 범위 검사는 유지합니다. 공급원은 기본 1 ms 선형 ramp로
0→10 kPa입니다. 기본 2 µm의 dt는 0.04 µs이며, 5 µm 비교 조건에서는
dt=0.25 µs입니다. 동일 물리 시간을 비교하려면 스텝 수를 다르게 지정해야 합니다.

```bash
make -j2
./geometry_validate 5 results/short_geometry_dx05
./geometry_validate 2 results/short_geometry_dx02
python3 preview_geometry.py

./main --outer-boundary open --reservoir-pressure-pa 10000 \
  --pressure-ramp-s 0.001 --steps 30000 \
  --diagnostic-every 1 --output-every 2000 \
  --output results/short_10kPa_ramp1ms_30000 \
  > results/short_10kPa_ramp1ms_30000.log 2>&1
python3 verify_vtk.py results/short_10kPa_ramp1ms_30000
python3 analyze.py
```

과거 결과는 이전 5 mm 영역의 계산이며 새 형상의 결과와 구분해야 합니다.
이전 README는 results/before_shortening_README.md에 보관했습니다.
축소 형상에서도 10 kPa 장시간 안정성을 보장하지 않습니다. 과거 1 ms ramp
실행은 3396스텝에 발산했으며 영역 축소 후에는 별도 검증이 필요합니다.

빌드 후 이전 5 µm에서 20스텝 단기 계산과 VTK 검사는 PASS입니다.
2 µm에서는 초기화 및 초기 VTK 검사를 통과했습니다(시간 적분은 수행하지 않음).
실제 converter의 dx/dt는 각각 5 µm/0.25 µs, 2 µm/0.04 µs로 확인했습니다.
이 결과는 축소 형상의 적용 확인이며 고압 장시간 안정성 검증은 아닙니다.

현재 기본 dx=2 µm, dt=0.04 µs입니다. 30,000스텝은 1.2 ms이며,
1 ms ramp 완료 시점은 25,000스텝입니다. 물성과 1 ms ramp는 유지합니다.


## 기본 2 µm 실행 시간 측정

`results/short_dx02_benchmark100/runtime_benchmark.json`에 단일 solver 프로세스,
매 스텝 진단 조건의 100스텝 실측을 저장했습니다. 초기화 및 VTK를 제외한
10–90스텝 처리 속도는 약 0.22577초/스텝(4.429스텝/초)이었습니다.
100스텝 전체 실행은 초기화와 두 VTK를 포함해 약 28.72초였습니다.

30,000스텝 순수 연산 예상은 약 6773초, **1시간 53분**입니다. 출력과 부하
변동을 감안하면 약 **2시간 전후(여유 있게 2–2.5시간)**를 예상할 수 있습니다.
이는 현재 장비의 짧은 실측 외삽이며 완료 시간이나 장시간 안정성을 보장하지
않습니다. 초기 100스텝 및 VTK 필드 검사는 PASS입니다.

현재 dx=2 µm에서 30,000스텝은 1.2 ms이며, 1 ms ramp는 25,000스텝에
10 kPa에 도달합니다. 5 µm일 때의 30,000스텝(7.5 ms)과 물리 시간이 다릅니다.

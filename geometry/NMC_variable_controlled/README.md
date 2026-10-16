# NMC_variable_controlled

`Cell_horizontal_wide_pore`에서 30,000 step을 통과한 실제 공기 기준 모델을 복사한 변수 제어용 폴더입니다. 원형 고체 배치와 계산 알고리즘은 유지했습니다. 기존 폴더의 결과는 복사하거나 덮어쓰지 않습니다.

## 실행

```bash
cd /home/soooo/OpenLB/release/openlb_electrolytefilling/geometry/NMC_variable_controlled
make -j2
./run.sh --electrolyte-pressure-pa 2500 --steps 30000
```

압력과 계산 스텝 수를 다음처럼 바꿉니다.

```bash
./run.sh --electrolyte-pressure-pa 5000 --steps 10000
./run.sh --electrolyte-pressure-pa 1000 --steps 20000 --output results/p1000_20k
./run.sh --help
```

압력 단위는 Pa이며, 상단·좌우 개방 경계의 0 Pa에 대한 공급원 압력입니다. 하단은 공급원 이외에 벽입니다. 지정 압력까지 기본 0.001초 동안 선형 증가합니다. 즉시 해당 압력을 가하는 조건은 `--pressure-ramp-s 0`입니다.

옵션을 생략하면 압력 2500 Pa, 30000 step입니다. `--reservoir-pressure-pa`도 같은 압력 옵션입니다. 같은 변수를 중복 지정하면 마지막 값이 적용됩니다. 스크립트는 옵션을 실행 파일에 전달하고 매 실행마다 `results/run_날짜_시간_PID`를 만듭니다. `--output`을 지정하면 그 경로를 사용하며, 기존 경로가 있으면 덮어쓰지 않고 종료합니다. 상대 경로는 이 폴더 기준입니다.

## 계산 시간의 의미

`--steps`는 총 계산 횟수이며 물리 시간 간격 Δt 자체가 아닙니다. 기본 조건에서 Δt = 2.5×10⁻⁷ s이므로:

| 스텝 수 | 물리 시간 |
|---:|---:|
| 2000 | 0.0005 s (0.5 ms) |
| 10000 | 0.0025 s (2.5 ms) |
| 30000 | 0.0075 s (7.5 ms) |

Δt = (τ_liquid − 0.5) Δx² / (3 ν_liquid)로 결정됩니다. `--dx-um`, `--tau`, `--nu-liquid` 변경 시 Δt도 변하며, phase mobility와 안정성도 달라질 수 있습니다. 따라서 스텝 수 변경과 물리 Δt 변경은 구별해야 합니다. 실제 적용값은 각 실행의 `run_config.csv`에 기록됩니다. 벽시계 연산 시간은 장비와 출력 주기에 따라 달라집니다.

## 유지된 기본 조건

- 격자 5 µm, 액체/상장 τ = 0.53, 접촉각 30°.
- 액체 밀도 1000 kg/m³, 동점성 1e-6 m²/s.
- 기체 밀도 1.2 kg/m³, 동점성 1.5e-5 m²/s, 온도 변화 미고려.
- 표면장력 1e-5 N/m: 기존 수치 검증값이며 실제 물–공기 값은 아닙니다.
- corrected RLB, bounded 물성 보간, conservative phi projection 유지.
- 너비 4.3 mm, 공급원 높이 0.5 mm, 그 위 다공성 높이 약 1.333 mm.
- 고정 원형 고체 직경 50–200 µm, 목표 공극률 0.40, 최소 입자 간격 20 µm, 동일 난수 seed.

phi projection을 적용한 안정화 모델입니다. projection 이전 phi의 범위 초과는 발생할 수 있으므로, 종료 없이 실행되었다는 사실만으로 물리 정확도나 국소 계면 거동이 검증된 것은 아닙니다. 압력 등의 변수를 바꾸면 기존 30,000 step 통과 결과가 보장되지는 않습니다.

기준 검증 기록: [기존 30,000 step 보고서](../Cell_horizontal_wide_pore/results/air_30k_v1/REPORT.md). 고체 관련 설명: [SOLID_PROPERTIES.md](SOLID_PROPERTIES.md).

## 출력과 중단

VTK는 기본 100 step마다, 진단 CSV는 매 step 저장됩니다. `--output-every 500`, `--diagnostic-every 10`으로 간격을 변경할 수 있습니다. 실행 로그는 출력 폴더 옆의 `.log`, VTK 검사 로그는 `.vtk_check.log`입니다. 정상 종료 후 스크립트가 `verify_vtk.py`를 실행합니다.

Ctrl+C 또는 오류로 중단되면 이미 디스크에 저장된 출력까지만 남습니다. 체크포인트 재시작 기능은 없으며 다음 실행은 초기 상태부터 시작합니다. 실행 중 진행 상황은 다른 터미널에서 `tail -f results/<실행이름>.log`로 볼 수 있습니다.

## 이 복사본의 확인

2500 Pa 및 5000 Pa 목표 압력으로 각각 10 step 실행하여 CLI 전달, 설정 기록, VTK 검사를 통과했습니다. 기본 조건의 첫 10 step 진단은 기존 30,000 step 실행과 완전히 일치합니다. 이 짧은 검사는 ramp 초기 구간의 동작 확인이며, 목표 압력 도달 후 안정성 검사는 아닙니다. 상세 결과는 `CLI_VALIDATION.json`에 기록했습니다. 새 폴더에서 장기 30,000 step 검증을 반복한 것은 아니며, 위 기준 보고서의 검증된 수치 계산 코드를 유지했습니다.

## 중력 추가

중력은 기본 0이며 `--gravity-y`로 부호를 포함한 y 가속도(m/s²)를 지정합니다. −y 방향의 지구 중력, 공급 압력 2500 Pa, 30000 step:

```bash
./run.sh --gravity-y -9.81 --electrolyte-pressure-pa 2500 --steps 30000
```

액체와 기체 유체 노드에 동일한 가속도를 적용하고 밀도에 따른 체적력이 작용합니다. 격자 가속도는 g Δt²/Δx로 환산하며 기본 격자에서 −1.22625e-7입니다. 공급원은 기존 압력 유지 방식이며 외곽은 기존의 0 Pa 개방 경계입니다. 중력은 시작부터 일정하게 적용하고 압력만 기존 ramp를 적용합니다. `run_config.csv`에 물리·격자 가속도를 기록합니다. 중력 조건의 장기 안정성은 기존 무중력 검증과 별도로 확인해야 합니다.

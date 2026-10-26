# NMC_gravity_duty_5to1

기존 연속 가진 1g 모델을 복사해, 주파수 10→100→10Hz의 부드러운 왕복 변화와 가진 5ms + 휴진 1ms를 적용한 모델이다.
기존 모델과 완료 결과는 수정하지 않는다.

| 항목 | 기본 조건 |
|---|---|
| 가진 진폭 | 1g = 9.81m/s² 편진폭 |
| 순간 주파수 | 10→100→10Hz, 코사인 형태로 부드럽게 왕복 |
| 주파수 왕복 시간 | 120ms, t=0/60/120ms에서 10/100/10Hz |
| 가진 / 휴진 | 5ms / 1ms, 5:1 |
| 세트 주기 / 듀티비 | 6ms / 5/6 (83.333%) |
| 총 계산시간 / 세트 수 | 120ms / 20세트 |
| 가진 누적시간 / 휴진 누적시간 | 100ms / 20ms |
| 휴진 중 조건 | 가진 가속도 0, 중력 −9.81m/s² 및 공급압력 2500Pa 유지 |
| 위상 | 휴진 중에도 전역 시간에 따라 연속 진행 |
| 구조 / 물성 / 경계 | 기존 모델과 동일 |
| 격자 / 시간 간격 | 5µm / 62.5ns |
| 총 스텝 | 1,920,000 |
| VTK / 진단 CSV | 10,000스텝마다 / 매 스텝 |

## 가진 파형

`t`를 초 단위 시간, `r = t mod 0.006`으로 두면 계산 구간 0 ≤ t < 0.12에서:

주파수 왕복 시간 `T=0.12초`일 때:

`f(t) = 55 − 45 cos(2πt/T)` [Hz]

`φ(t) = 2π × 55 × t − 45 × T × sin(2πt/T)` [rad]

위상은 순간 주파수의 시간 적분으로 계산한다. `2πf(t)t`로 대체하지 않는다.

- 0 ≤ r < 0.005: `a_excitation_y = M × 9.81 × sin(φ(t))`
- 0.005 ≤ r < 0.006: `a_excitation_y = 0`
- `a_effective_y = −9.81 + a_excitation_y`

구간은 [시작, 끝)이며, 최종 t=0.12초에서 가진은 0이다.
첫 가진은 0–5ms, 첫 휴진은 5–6ms, 두 번째 가진은 6–11ms,
두 번째 휴진은 11–12ms이다. 마지막은 가진 114–119ms, 휴진 119–120ms이다.

주파수와 위상은 휴진 중에도 전역 시간에 따라 계속 진행하며 세트 시작 시 초기화하지 않는다.
고정 주파수 정현파와 달리, 가속도 파형의 순간 주파수가 연속적으로 변한다.
주파수는 0ms:10Hz, 30ms:55Hz, 60ms:100Hz, 90ms:55Hz, 120ms:10Hz이다.
120ms 동안 주파수는 한 번 왕복하고, 가진·휴진은 20세트 반복한다.
전체 반송파의 위상 진행은 6.6주기다. 실제 가진 중 주기 상당량은 휴진을 제외해 따로 적분한다.
가진 창의 길이가 5ms로 고정되어 있으므로 각 창이 포함하는 사인파 주기의 비율은 시간에 따라 다르다.

휴진 전환에는 별도 진폭 램프를 적용하지 않아 가속도가 불연속일 수 있다.
10~100Hz는 반송파의 **순간 주파수 범위**다. 가진·휴진 게이트를 적용한 최종 파형의
푸리에 스펙트럼이 이 대역에만 제한된다는 뜻은 아니다.

고체 구조는 고정이고, 유체에 외력 가속도를 적용한다.
휴진은 계산 정지가 아니며 유동 계산은 계속 진행한다.

## 빌드와 실행

OpenLB가 설치된 `release` 폴더 안에 이 저장소를 둔다. 저장소 루트는
`release/openlb_electrolytefilling`이며, `geometry` 폴더가 바로 아래에 있다.
기존에 OpenLB와 MPI 실행 환경을 준비한 조원은 `release` 폴더에서 다음을 실행한다.

```bash
cd openlb_electrolytefilling
git pull origin main
cd geometry/NMC_gravity_duty_5to1
make -B -j4
./run_gravity_scale.sh 1
```

저장소를 처음 받는 경우 OpenLB `release` 폴더에서:

```bash
git clone https://github.com/sooyeonji99/openlb_electrolytefilling.git
cd openlb_electrolytefilling/geometry/NMC_gravity_duty_5to1
make -B -j4
./run_gravity_scale.sh 1
```

필요 환경은 C++20 지원 컴파일러(GCC 12 이상), GNU Make, OpenMPI의 `mpic++`·`mpirun`, Python 3이다.
OpenLB 소스는 이 저장소에 포함하지 않는다. 이번 실행에 사용한 OpenLB revision은
`882924a9cfc8dcdf82909e39530789cb6f5ebcc4` (`1.9.0-21-g882924a9`)이다.
OpenLB 루트의 `external/lib/libz.a`와 `libtinyxml2.a`가 필요하다.
없다면 OpenLB `release` 폴더에서 `make -C external lib`,
`make -C external zlib tinyxml2`를 순서대로 실행한다.
별도 위치에 설치된 OpenLB를 사용할 때는 `make -B -j4 OLB_ROOT=/path/to/OpenLB`로 지정한다.

첫 인자는 g 배율, 둘째는 선택적 출력 경로, 셋째는 MPI 프로세스 수(기본 4)다.
기본 출력은 `results/gravity_sweep_10to100_duty_5to1_1g_날짜_시간`이며 기존 결과를 덮어쓰지 않는다.
이 명령은 120ms 전체 계산을 시작한다. 2026-10-26에 MPI 4개 프로세스로 전체 계산을 시작했다. 시작 구간의 수치 유효성과 이전 모델과의 형상 일치를 확인했으며, 전체 계산 완료 및 장시간 안정성은 아직 확인하지 않았다.

직접 옵션 지정:

```bash
./run_mpi.sh --ranks 4 \
  --outer-boundary top-open --dx-um 5 --duration-s 0.12 \
  --electrolyte-pressure-pa 2500 --pressure-ramp-s 0 \
  --gravity-y -9.81 --gravity-scale 1 \
  --excitation-min-frequency-hz 10 --excitation-max-frequency-hz 100 \
  --excitation-sweep-period-s 0.12 \
  --excitation-on-s 0.005 --excitation-rest-s 0.001 \
  --tau 0.5075 --phase-tau 0.5075 \
  --output-every 10000 --diagnostic-every 1 \
  --output results/duty_5to1_1g
```

`--excitation-on-s`와 `--excitation-rest-s`의 단위는 초다. 휴진 시간을 0으로 설정하면 계산 구간 내 연속 가진이 된다.
Makefile은 이 모델 전용 MPI 코어를 `.build/mpi`에 빌드한다.
다른 PC로 옮기면 `make -B -j4`로 다시 빌드한다.

## 결과 및 검증

- `run_config.csv`: 주파수 최솟값·최댓값·왕복 시간, 가진·휴진 시간, 듀티비, 세트 수, 총 가진시간, 위상 정책.
- `forcing_history.csv`: 매 스텝 가진 상태 `excitation_active`, 가진 가속도, 중력 포함 총 가속도, 압력, 순간 주파수, 적분한 반송파 위상.
- `analysis_summary.json`: 계산 완료 및 수치 유효성, 가진 누적시간, 반송파 위상 주기 수.
- `results/forcing_schedule_1g.csv`: 전체 120ms 처방 파형. 유동 계산 결과는 아니다.
- [VALIDATION.md](VALIDATION.md): 이번 검증의 범위와 결과.
- `results`와 실행파일은 Git에서 제외하며 각 PC에서 새로 생성한다.

작성자 PC의 기존 연속 가진 데이터 보관 위치(Git에 포함하지 않음):
`../../archives/NMC_gravity_scale_continuous_1g_20261022_185042/`
보관본은 결과 전체, 실행 로그, 실행파일, 모델 소스·설정을 포함하며
`ARCHIVE_MANIFEST.json`에 모든 파일의 SHA-256 대조 결과가 기록된다.

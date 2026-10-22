# NMC_forcesweep_longterm

하단 중앙의 electrolyte 공급원에서 상부 개방면으로 침투·배출되는 2차원 다상 LBM 모델이다. 기존 계산 폴더 및 결과와 독립적이다.

## 구조와 물성

|항목|설정|
|---|---|
|가로|기존과 같은4.3 mm|
|전체 높이|3.670 mm: 기존 격자상 높이1.835 mm의2배|
|공급원|하단 중앙1 mm ×0.5 mm, 공급원 높이는 유지|
|외곽|좌우·하단 벽, 상부0 Pa 게이지압 개방. 상부로 액체와 기체 모두 배출 가능|
|공급 압력|시작부터 종료까지2500 Pa 게이지압 유지. 초기ramp와 중간 압력 해제 없음|
|기공률|목표0.4. 실제 격자값은geometry_validation.csv에 기록|
|입자|원형50–200 µm, 최소 간격20 µm, seed20261014. 확장 영역에 재배치|
|격자|5 µm,861×735=632,835개 격자점|
|시간 간격|6.25e-8초, 액체tau=phase tau=0.5075|
|액체|rho1000 kg/m³, nu1e-6 m²/s|
|기체|rho1.2 kg/m³, nu1.5e-5 m²/s|
|표면장력·접촉각|1e-5 N/m,30도|
|중력|전 구간y방향−9.81 m/s²|
|가진 진폭|기존의 전체 높이2배 규칙 유지. 변위 편진폭A=7.34 mm|
|출력|VTK10,000step마다=0.625ms, 진단CSV매step|

## 가진 일정

|물리시간|조건|
|---|---|
|0–0.030초|무가진|
|0.030–0.530초|10→100 Hz 선형 증가, +180 Hz/s|
|0.530–1.030초|100→10 Hz 선형 감소, −180 Hz/s|
|1.030–1.060초|무가진|

무가진 중에도 중력과2500Pa 공급은 유지된다. 총1.06초=16,960,000step. 스윕은 총55주기(증가27.5주기, 감소27.5주기)이다.

순간 주파수f(t)를 적분한 연속 위상 `θ(t)=2π∫f(t)dt`를 사용한다. `sin(2π f(t)t)`로 대신하지 않는다. 변위 `x=A sinθ`의2차 미분에서 구한 관성 가속도는:

`a_inertial = A[(2πf)^2 sinθ − 2π(df/dt) cosθ]`

이 값에 중력을 더해 고정된 다공체 좌표계의 외력 가속도로 적용한다. 실제 고체 이동·변형은 계산하지 않는다. 시작·종료 시 속도 불연속에 의한 충격은 포함하지 않으며, 주파수의 선형 증가/감소가 만나는 지점에서는df/dt의 부호 변화로 가속도에 단차가 생긴다. 휴진 중 유체 운동을 강제로 정지시키지 않는다.

## 실행

```bash
cd /home/soooo/OpenLB/release/openlb_electrolytefilling/geometry/NMC_forcesweep_longterm
./run_forcesweep.sh 4
```

6개 프로세스는 `./run_forcesweep.sh 6`으로 실행한다. 출력 폴더를 직접 지정하려면:

```bash
./run_forcesweep.sh 4 results/forcesweep_1060ms
```

기본 출력은 날짜·시각을 붙인 새 폴더이며 기존 결과를 덮어쓰지 않는다. ParaView에서는VTM/PVD를 연다. 다른PC에서는OpenLB설치 경로를Makefile의OLB_ROOT에 맞춘 후 `make -j2`로 빌드한다.

직접 옵션을 바꾸는 예:

```bash
./run_mpi.sh --ranks 4 \
  --outer-boundary top-open --duration-s 1.06 \
  --electrolyte-pressure-pa 2500 --pressure-ramp-s 0 --gravity-y -9.81 \
  --pre-rest-s 0.03 --sweep-up-s 0.5 --sweep-down-s 0.5 --post-rest-s 0.03 \
  --sweep-low-hz 10 --sweep-high-hz 100 --excitation-height-factor 2 \
  --tau 0.5075 --phase-tau 0.5075 \
  --output-every 10000 --diagnostic-every 1 \
  --output results/forcesweep_1060ms
```

짧은 검증은 `--duration-s` 대신 `--steps N`으로 제한할 수 있다. 두 옵션은 함께 지정하지 않는다. `--duration-s`를 생략하면 지정한4개 구간 길이의 합을 계산한다. 주파수와 구간 길이 옵션은위 명령 참조.

## 기록과 검증

- `forcing_history.csv`: 압력·변위·속도·가속도와 함께 순간 주파수·주파수 변화율·위상·구간명을 저장한다.
- `multiphase_diagnostics.csv`: 시간별 기공 평균phi, 범위, 보정량, 수치 유효성 등.
- `outlet_diagnostics.csv`: 상부 출구 속도·역류·비유한 값 등.
- `results/forcing_schedule_full.csv`: 전체1.06초의 처방 파형. 유동 계산 결과가 아니다.
- `results/VALIDATION_REPORT.md`: 단시간 검증 결과와 범위.

공급원phi를1로 유지하므로 상부에서 액체가 배출돼도 계속 공급된다. 기체도상부로 배출되며 역류 시 외부 기체를 적용한다. 단시간검증 PASS와 전체1.06초 안정성은 별개다. 액체 점성은기존 물 수준값이고 표면장력은기존 수치 검증값이므로 실제 물/전해액 물성 전체를 재현하는 설정은 아니다. 장시간실행은 자동시작하지 않는다.

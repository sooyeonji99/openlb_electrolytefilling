# NMC_gravity_scale

100 Hz의 y방향 사인파 가진을 처음부터 가하며, 가진 가속도 진폭을 중력 단위 배율로 조절하는 다상 LBM 모델이다. 기존 폴더와 결과는 변경하지 않는다.

## 조건

|항목|설정|
|---|---|
|전체 크기|가로4.3 mm × 세로3.67 mm|
|격자 간격|5 µm,861×735개 격자점|
|기공 구조|목표 기공률0.4, 원형 입자50–200 µm, 최소 간격20 µm, seed20261014|
|공급원|하단 중앙1 mm ×0.5 mm, 액체상 유지|
|공급 압력|시작부터 종료까지2500 Pa 게이지압, ramp없음|
|외곽|좌우·바닥 벽, 윗면0 Pa 게이지압 개방|
|상부 배출|액체·기체 배출을 허용. 역류 시 외부 기체 적용|
|중력|y방향−9.81 m/s²|
|가진|t=0부터100 Hz 연속, 위상0에서 시작|
|시간|120 ms=0.12초,100 Hz12주기|
|시간 간격|6.25e-8초, 목표1,920,000step|
|액체|rho1000 kg/m³, nu1e-6 m²/s|
|기체|rho1.2 kg/m³, nu1.5e-5 m²/s|
|액체·상장tau|0.5075·0.5075|
|접촉각·표면장력|30도·1e-5 N/m, 기존 수치 검증값 유지|
|출력|VTK10,000step마다(0.625ms), 진단CSV매step|

고체는 고정된 원형 장애물이다. 실제 고체 이동·변형을 계산하는 대신, 가속도를 유체의 외력 항으로 적용한다.

## 가진 magnitude

`--gravity-scale M`으로 가속도 편진폭을 정한다.

`a_excitation(t) = M × 9.81 × sin(2π × 100 × t)`

`a_effective_y(t) = gravity_y + a_excitation(t)`

|M|가진 가속도 범위|중력−9.81을 포함한 총 가속도 범위|
|---|---|---|
|0|0|−9.81 m/s²|
|0.5|±4.905|−14.715~−4.905 m/s²|
|1|±9.81|−19.62~0 m/s²|
|2|±19.62|−29.43~9.81 m/s²|

배율의 기준g는9.81 m/s²이다. `--gravity-y`는 정적 중력 항만 바꾸며 이 기준g와 독립적이다. 변위 진폭을 고정하는 설정이 아니며 변위·속도 경로는 지정하지 않는다. 휴진이나 주파수 스윕은 없다. M=0은 동일 압력·중력·구조의 무가진 기준 케이스다.

## 4코어 실행

```bash
cd /home/soooo/OpenLB/release/openlb_electrolytefilling/geometry/NMC_gravity_scale

# 기본1g,4프로세스 MPI,120ms
./run_gravity_scale.sh 1

# 2g,출력 폴더 지정
./run_gravity_scale.sh 2 results/gravity_2g

# 0.5g
./run_gravity_scale.sh 0.5 results/gravity_0p5g

# 무가진 기준
./run_gravity_scale.sh 0 results/gravity_0g
```

첫 인자=magnitude배율, 두 번째=출력 경로(생략하면 날짜·시각 포함 새 폴더), 세 번째=MPI프로세스 수(기본4)다. 기존 결과를 덮어쓰지 않는다.

직접 옵션을 지정하는 명령:

```bash
./run_mpi.sh --ranks 4 \
  --outer-boundary top-open --duration-s 0.12 \
  --electrolyte-pressure-pa 2500 --pressure-ramp-s 0 \
  --gravity-y -9.81 \
  --excitation-frequency-hz 100 --gravity-scale 1 \
  --tau 0.5075 --phase-tau 0.5075 \
  --output-every 10000 --diagnostic-every 1 \
  --output results/gravity_1g
```

`--gravity-scale`만 바꾸면 magnitude를 비교할 수 있다. 격자·시간 간격·물성은 함께 유지한다. `--duration-s` 대신 `--steps N`으로 짧게 계산할 수 있으며 두 옵션은 동시에 지정하지 않는다.

다른PC에서는Makefile의OLB_ROOT를설치경로에 맞추고 `make -j4`로 다시 빌드한다. CPU에 맞춘 빌드이므로 다른PC로 실행파일만 복사하지 않는다. 실행이전 `mpirun`과OpenLB의MPI 빌드가 필요하다.

## 결과와 검증

- `multiphase_diagnostics.csv`: 충전율·phi범위·속도·유한성 등.
- `forcing_history.csv`: 공급압력·가진 가속도·중력 포함 총 가속도·주파수·위상·배율.
- `outlet_diagnostics.csv`: 상부 출구·역류·비유한 값 진단.
- `multiphase_validation.csv`, `analysis_summary.json`: 종료 판정과 원인.
- `vtkData/cell_horizontal.pvd`: ParaView에서 여는 전체 시계열. VTM/VTI는 함께 보관한다.
- `results/forcing_schedule_1g.csv`: 120ms의 처방 파형으로, 유동 계산 결과가 아니다.
- `results/VALIDATION_REPORT.md`: 짧은 계산 검증과 한계.

장시간120ms 계산은 자동 시작하지 않는다. 짧은 검증의PASS가 큰M에서의 장시간 안정성이나120ms 내 상부 액체 도달을 보장하지는 않는다.

# NMC_timecost

8.6 mm 직경 원의 1/4 부채꼴 NMC 단면에 대한 형상 및 pore-resolved LBM
계산비 검증 모델입니다. 기존 NMC 및 `Example_validation` 파일은 수정하지
않습니다.

현재 형상 해석은 다음과 같습니다.

- 외곽 원 반지름: 4.3 mm
- 중심 electrolyte 원의 **직경**: 1.785 mm
- 따라서 중심 electrolyte 반지름: 0.8925 mm
- 중심 영역 바깥: solid와 pore가 공존하는 annulus
- 목표 공극률: 기존 NMC 설정과 같은 0.30

입도분포가 아직 지정되지 않았으므로 외곽 다공성 영역은 50 µm pitch의 규칙적
원형 solid 배열로 생성합니다. 반지름은 목표 공극률에서 계산되며 약 23.60 µm
입니다. 이 배열은 형상 처리, 연결성 및 계산비 검증용입니다. 실제 NMC 입자
분포를 의미하지 않습니다.

## 실행

```bash
cd /home/soooo/OpenLB/release/openlb_electrolytefilling/geometry/NMC_timecost
make
./nmc_timecost
make multiphase_smoke -j2
./multiphase_smoke
python3 analyze.py
```

별도 결과 폴더와 격자를 지정할 수 있습니다.

```bash
./nmc_timecost --dx-um 2 --porosity 0.30 \
  --outer-diameter-mm 8.6 --electrolyte-diameter-mm 1.785 \
  --particle-pitch-um 50 --time-s 10 --output results/run_02
```

재실행하면 지정한 결과 폴더의 파일을 갱신합니다.

## 출력

- `quarter_sector_geometry.vtk`: ParaView 형상
- `geometry_validation.csv`: 면적, 공극률, 연결성
- `validation_result.csv`: 자동 PASS/FAIL
- `time_cost.csv`: 해상도별 셀, Δt, 스텝, 메모리 및 런타임 추정
- `results/multiphase_smoke/multiphase_diagnostics.csv`: OpenLB 시간 이력
- `results/multiphase_smoke/multiphase_validation.csv`: multiphase PASS/FAIL

VTK material 번호는 0=부채꼴 외부, 1=초기 electrolyte, 2=pore,
3=solid입니다.

## 현재 검증 결과

기본 `dx=2 µm` 형상은 다음 결과로 **PASS**입니다.

- 격자 크기: 2151 × 2151
- 부채꼴 내부 셀: 3,632,607
- 외곽 면적 상대오차: 0.0580%
- 중심 electrolyte 면적 상대오차: 0.2868%
- annulus 공극률: 0.296286
- 중심부와 연결된 pore 비율: 0.999967
- pore가 외곽 원호까지 연결됨

## 시간비용 해석

완화시간 0.8, 동점도 `1e-6 m2/s`이면
`dt=(tau-0.5) dx^2/(3 nu)`입니다. `dx=2 µm`에서는 `dt=4e-7 s`이고
10초에 25,000,000 스텝이 필요합니다. 실제 `dx=5 µm` smoke run에서 측정한
격자 노드당 메모리를 위로 반올림한 500 byte/node를 적용하면 `dx=2 µm`의
정사각형 OpenLB 할당 영역은 약 2.15 GiB입니다. 50 MLUPS에서도 약 504시간이
필요합니다. 실제 curved wall,
통신, VTK 출력과 진단 비용은 이 추정보다 늘어납니다.

기존 최소 입자 반경 0.1 µm까지 직접 해상하려고 `dx=0.2 µm`를 사용하면 약
3.63억 부채꼴 셀, 약 215 GiB의 OpenLB 할당 메모리, 25억 스텝이 필요합니다.
현재 장비의 15 GiB RAM에서는
실행할 수 없습니다. 따라서 전체 8.6 mm 형상에는 다음 중 하나가 필요합니다.

1. 공극 평균화 Darcy/상대투과도 모델
2. 대표 체적요소(RVE)에서 pore-scale 물성을 보정한 다중스케일 모델
3. 영역 분할 MPI와 대용량 메모리를 사용하는 HPC 계산

## OpenLB multiphase smoke 검증

`multiphase_smoke.cpp`는 앞선 validation에서 통과한 conservative phase-field,
Liang coupling, incompressible BGK 조합을 같은 8.6 mm 부채꼴에 적용합니다.
중심 반지름 0.8925 mm 내부는 `phi=1`, 외곽 pore는 `phi=0`, solid와 외곽은
60도 wetting wall입니다. 계산비를 제한하기 위해 smoke run은 `dx=5 µm`,
표면장력 `1e-5 N/m`, 200스텝만 사용합니다.

- OpenLB `SuperGeometry`: 정상
- 물리시간: 0.0005 s
- 최종 phi 범위: `3.33e-16`–`1.00002`
- 최대 격자속도: `0.000609472`
- 상태: **PASS**
- 측정 wall time: 9.97초
- 측정 최대 RSS: 355,468 KiB

이는 전체 장시간 침투 결과가 아니라 재료 매핑, 초기 계면, wetting wall,
NS/phase-field coupling과 출력이 실제 OpenLB에서 작동한다는 통합 smoke
검증입니다. 정량적인 실제 NMC 침투 계산은 입도분포와 공급 경계가 확정된 뒤
RVE 또는 공극 평균화 전략으로 구성해야 합니다.

## 권장 계산 구성

전체 8.6 mm 부채꼴은 공극 평균화 Darcy 모델로 계산하고, 100–300 µm 크기의
대표 체적요소(RVE)에만 현재 multiphase LBM을 적용합니다. RVE에서 투과도,
모세관압 및 상대투과도를 구해 전체 모델에 전달하면 형상을 유지하면서 계산
시간을 크게 줄일 수 있습니다. 10 Hz 가진은 전체 모델에서 1초(10주기)를 먼저
계산하고, 계면 진동의 주기 평균값과 누적 침투량을 함께 비교합니다.

## MPI 병렬 계산

현재 `/home/soooo/OpenLB/olb-release/config.mk`는 `PARALLEL_MODE := OFF`인
직렬 빌드입니다. MPI 실행 전 OpenLB 루트에서 MPI 설정으로 라이브러리를 다시
빌드해야 합니다.

```bash
cd /home/soooo/OpenLB/olb-release
cp config.mk config.serial.mk
cp config/cpu_gcc_openmpi.mk config.mk
make clean
make -j12
```

계산 코드의 cuboid 수는 MPI rank 수와 같아야 합니다. `multiphase_smoke.cpp`는
이미 `singleton::mpi().getSize()`를 사용하도록 구성되어 있어 직렬에서는 1개,
MPI에서는 rank 수만큼 분할합니다. 다만 현재 이 디렉터리의 간이 `Makefile`에
있는 `multiphase_smoke` target은 직렬 전용이므로, 본 계산 solver는 OpenLB의
표준 `default.mk` 방식으로 MPI 설정을 상속해 빌드해야 합니다. 빌드 후에는
다음처럼 실행합니다.

```bash
mpirun -np 4 ./solver
```

먼저 `-np 1`, `-np 4`, `-np 8`로 같은 짧은 계산을 수행해 결과가 같은지와
wall time을 비교합니다. 출력 주기를 길게 두면 MPI 파일 입출력 병목을 줄일 수
있습니다. 작업이 끝난 뒤 직렬 설정은 `cp config.serial.mk config.mk` 후
`make clean && make -j12`로 복원할 수 있습니다.

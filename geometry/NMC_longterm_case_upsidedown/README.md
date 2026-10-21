# NMC_longterm_case_upsidedown

기존 NMC_longterm_case의 구조를 상하 반전한 2차원 pore-resolved 다상 LBM 모델이다. 기존 실행 및 결과와 독립된 폴더이며 기존 세 케이스 순차 계산은 계속 진행한다.

- 공급원: 윗벽 중앙, 가로1 mm × 세로0.5 mm. 액체상 유지, 게이지압0 Pa.
- 경계: 위 벽, 좌우 및 아래0 Pa 개방. 고체/벽에 기존 접촉각30도 적용.
- 기존 seed20261014의330개 원형 입자를 동일하게 상하 반전한다. 입자 지름50–200 µm, 최소 간격20 µm. 재무작위 생성하지 않는다.
- 목표 porosity0.4; 실제 격자 기공률은 약0.398094이며 geometry_validation.csv에서 확인한다.
- 크기:4.3 mm × 격자상1.835 mm. dx5 µm,861×368격자.
- 중력: lab 좌표의 y방향−9.81 m/s² 유지. 공급원이 위쪽으로 옮겨져 공급원에서 기공으로의 진행 방향과 중력의 관계가 기존과 반대가 된다.
- 추가 압력 가압 없이100 Hz 연속 가진을 t=0부터 적용. A=모델 높이2배=3.67 mm, 관성 가속도 진폭 약1448.86 m/s². 중력과 젖음에 따른 자연침투도 함께 발생하며 이를 물리적으로 제거하지 않는다.
- dt6.25e-8초, liquid/phase tau0.5075, liquid rho1000·nu1e-6, gas rho1.2·nu1.5e-5, sigma1e-5 N/m 유지.
- 목표60 ms=960,000 step=100 Hz6주기. VTK10,000 step마다(0.625 ms), 진단 매step.

## 실행

```bash
cd /home/soooo/OpenLB/release/openlb_electrolytefilling/geometry/NMC_longterm_case_upsidedown
./run_100hz_60ms.sh 4
```

첫 인자는MPI프로세스 수(기본4). 기존 계산과 동시에 실행하면 CPU를 나눠 쓰므로 속도가 느려질 수 있다. 이 환경 준비 과정에서는 짧은 검증만 실행하며60ms 장기 계산은 자동 시작하지 않는다.

출력 폴더를 지정하려면:

```bash
./run_100hz_60ms.sh 4 results/upsidedown_100hz_60ms
```

직접 옵션을 조절할 때:

```bash
./run_mpi.sh --ranks 4 \
  --outer-boundary open --steps 960000 \
  --electrolyte-pressure-pa 0 --gravity-y -9.81 \
  --excitation-frequency-hz 100 --excitation-height-factor 2 \
  --excitation-start-s 0 --excitation-duration-s 0.06 \
  --excitation-on-s 0 --excitation-off-s 0 \
  --tau 0.5075 --phase-tau 0.5075 \
  --output-every 10000 --diagnostic-every 1 \
  --output results/upsidedown_100hz_60ms
```

기존 출력은 덮어쓰지 않는다. 실패 시에도 로그·이상VTK·분석을 저장한다. ParaView에서는 VTM/PVD를 연다. 공급원은 상단에 붙은 유한 높이의 액체 저장 영역이다. 압력0 Pa에서도 액체상을 유지하므로 닫힌 유한 액적과는 다르다.

구조의 상하 반전은 전체 격자점에서 기존 모델과 대조한다. 짧은 검증의 PASS는60ms 안정성을 보장하지 않는다.

검증 완료: 전체316,848격자점 상하 반전 일치, 2프로세스 MPI100step, 공급원phi=1·0Pa·상단 벽/좌우하단 개방·VTK 검사 PASS. `results/VALIDATION.json` 참고.

# Cell_horizontal_wide_pore — 공기 물성 30000스텝 검증

**30000스텝(7.5 ms) 계산과 저장 VTK 301개 검사를 완료했습니다.**
결과: [검증 보고서](results/air_30k_v1/REPORT.md),
[요약 JSON](results/air_30k_v1/accepted_summary.json),
[실패 후보 비교](results/air_30k_v1/comparison.csv).

## 현재 기본 물성·조건

| 항목 | 값 |
|---|---|
| 기체 모델 | 상온·대기압 공기의 고정 참조 물성, 비압축성 근사 |
| 공기 밀도 | 1.2 kg/m³ |
| 공기 동점성 / 점성계수 | 1.5e-5 m²/s / 1.8e-5 Pa·s |
| 액체 밀도 / 동점성 | 1000 kg/m³ / 1e-6 m²/s |
| 표면장력 | 기존 검증용 1e-5 N/m |
| 명목 접촉각 | 30° |
| 공급 압력 | 0→2500 Pa, 1 ms 선형 ramp 후 유지 |
| 외곽 | 하단 벽, 상단·좌우 0 Pa 기준 압력 개방 |
| NS 액체 / 기체 tau | 0.53 / 0.95 |
| 상장 tau / 기본 mobility 계수 | 0.53 / 1e-6 m²/s |
| 격자 / 시간 간격 | 5 µm / 0.25 µs |
| 기본 실행 | 30000스텝, 진단 매 스텝, VTK 100스텝마다 |

온도/에너지 방정식과 압축성 공기 EOS는 추가하지 않았습니다.
공기 물성만 실제 대표값으로 변경했고, 액체–기체 표면장력을 실제 물–공기 값으로
교체하지 않았습니다. 고체는 [고정·비변형 벽 모델](SOLID_PROPERTIES.md)입니다.

형상은 4.3 × 1.835 mm, 861 × 368 노드, 공급원 1 × 0.5 mm,
직경 50–200 µm 원형 입자 330개, 최소 입자 간격 약 20 µm,
실제 공극률 0.398094입니다. 기공 연결 및 입자 형상 검사도 통과했습니다.

## 안정화를 위해 적용한 계산 방법

- NS: OpenLB 고밀도비 예제에 기반한 보정 정규화 충돌 모델 `crlb`.
  `physical_range_coupling.h`에 보정 밀도 기울기/Laplacian과 계면력 및 상장 flux 보정을 구현했습니다.
- 밀도·점성 환산: `--material-interpolation bounded`로 물리 범위 보호.
- 상장: `--phase-projection conservative`로 pore+출구의 총 phi를 보존하며
  diffuse interface에 보정량을 재분배합니다. 공급원/고체는 재분배 대상에서 제외합니다.
  분포함수 보정은 1차 모멘트를 유지하며, 기존 NS/상장 tau와 mobility 계수는 유지했습니다.
- 기존 phi/양의 밀도/NS tau/유한성 중단 검사를 유지했습니다. 속도 크기 중단은 기존처럼 비활성화입니다.

**이번 통과는 phi 범위 보정에 의존합니다.** 보정 전 phi는 순간적으로
-0.508~1.165까지 벗어났습니다. 보정 이후 범위는 반올림 오차를 제외하면 [0,1]이고,
최대 NS 격자 속도는 0.051545, 저장 VTK 최대값은 0.043269였습니다.
보정 연산의 최대 총량 잔차는 약 7.3e-9 phi×cell입니다.

이는 수치 진행 검증입니다. 원격 계면 사이 재분배가 가능하므로 국소 계면/기포 변화,
공급·배출을 포함한 질량 수지, 접촉각 및 모세관 압력의 정량 정확성을 별도로 검증해야 합니다.

## 실행

```bash
make -j2
./run.sh 30000
```

`run.sh`는 위 검증 조건을 명시적으로 사용하고 매번 새 결과 폴더를 만듭니다.
출력: `results/air_theta30_p2500_30000_<실행시각>_<PID>/` 및 같은 이름의 `.log`.
정상 종료 후 VTK 검사를 자동 실행합니다. 각 run_config.csv가 실제 실행 조건입니다.

개별 결과 검사:

```bash
python3 verify_vtk.py results/<실행폴더>
```

완료된 최초 공기 검증 보고서 재생성은 `python3 analyze_air_validation.py`입니다.
이 분석기는 `results/air_30k_v1/air_corrected_rlb`를 검사하며 최신 수동 실행을 자동 선택하지 않습니다.

## 과거 검증

- [0 Pa / 표면장력 대조군](results/zero_pressure_controls_v1/REPORT.md): 과거 접촉각 60°, 기체 밀도 100 kg/m³.
- [접촉각·기체 점성 대조군](results/wetting_gas_theta30_v1/REPORT.md): 과거 기체 밀도 100 kg/m³.
- 기존 validation_2000 및 최초 10000스텝 결과도 과거 물성/충돌 조건입니다.
- 과거 재현 스크립트는 해당 밀도·상장 보정 off·선형 물성 환산·BGK를 명시합니다.
- 이전 README와 실행 당시 소스는 results/air_30k_v1에 보관했습니다. 원본 OpenLB 라이브러리는 수정하지 않았습니다.

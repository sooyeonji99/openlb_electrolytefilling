# 완료된 계산 결과

두 실행 모두 동일 Geometry에서 총 1,920,000스텝(120ms)을 완료했으며 최종 검증은 PASS다.

| 실행 | 가진 조건 | 최종 기공 평균 액체상 비율 |
|---|---|---|
| Long term (비교 기준) | 100Hz 고정, 1g, 120ms 연속 가진 | 48.7129% |
| Sine sweep | 10→100→10Hz, 왕복 120ms, 1g, 가진 5ms·휴진 1ms | 48.7344% |

## 결과 위치

- Long term: [결과 요약](geometry/NMC_gravity_scale/results/gravity_scale_1g_20261022_185042/analysis_summary.json), [실행 조건](geometry/NMC_gravity_scale/results/gravity_scale_1g_20261022_185042/run_config.csv)
- Sine sweep: [결과 요약](geometry/NMC_gravity_duty_5to1/results/gravity_sweep_10to100_duty_5to1_1g_20261026_160905/analysis_summary.json), [실행 조건](geometry/NMC_gravity_duty_5to1/results/gravity_sweep_10to100_duty_5to1_1g_20261026_160905/run_config.csv)

각 폴더에 전체 시계열 CSV, 형상, 검증 결과, ParaView용 VTK 시계열을 포함한다.
상위 results 폴더의 실행 로그·분석 로그·병렬 설정도 포함한다.
실행파일 및 이전 데이터의 중복 보관본은 포함하지 않는다.
파일 수와 용량은 [RESULTS_MANIFEST.json](RESULTS_MANIFEST.json)에 기록한다.

## 대용량 결과 받기

두 결과는 합계 약 19.2GB이다. 대용량 CSV와 VTI는 Git LFS로 저장한다.
실제 업로드에는 저장소의 LFS 저장공간·전송량 한도가 충분해야 한다.

Git LFS를 설치한 PC에서 저장소 루트에서 실행한다.

```bash
git lfs install
git pull origin main
git lfs pull
```

처음에는 소스와 메타데이터만 받고 원본 데이터는 선택적으로 받으려면:

```bash
GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/sooyeonji99/openlb_electrolytefilling.git
cd openlb_electrolytefilling
git lfs install
# Sine sweep만 다운로드
git lfs pull --include="geometry/NMC_gravity_duty_5to1/results/gravity_sweep_10to100_duty_5to1_1g_20261026_160905/**"
# Long term만 다운로드
git lfs pull --include="geometry/NMC_gravity_scale/results/gravity_scale_1g_20261022_185042/**"
```

LFS 다운로드를 하지 않은 파일은 포인터 텍스트이므로 ParaView나 CSV 도구에서 직접 열 수 없다.
다운로드 후 각 결과 폴더의 `vtkData/cell_horizontal.pvd`를 ParaView에서 연다.

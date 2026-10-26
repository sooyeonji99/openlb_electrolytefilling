# 검증 기록 (2026-10-26)

## 처방 파형 및 빌드

- 전체 1,920,000스텝에서 가진/휴진 상태를 정수 스텝 기준과 대조: PASS.
- 가진 5ms + 휴진 1ms, 총 20세트; 가진 총 100ms, 휴진 총 20ms.
- 순간 주파수 10→100→10Hz, 왕복 시간 120ms.
- 주파수를 적분한 위상을 독립 중점 수치적분 및 위상 미분으로 검증: PASS.
- 12,001개 시점의 C++ 처방 파형과 Python 검증식 대조: PASS.
- 진폭 배율, 휴진 0, 잘못된 입력 검증: PASS.
- MPI 실행파일 빌드, Python·셸 문법 검증: PASS.

파형 검증 재현(모델 폴더에서 실행):

```bash
mkdir -p results
g++ -O2 -std=c++20 test_gravity_excitation.cpp -o test_gravity_excitation
./test_gravity_excitation
```

`results/forcing_schedule_1g.csv`는 처방 파형이며 유동 결과가 아니다.

## 기존 모델과의 동일성 및 실행 상태

2026-10-26 16:09 KST에 MPI 4개 프로세스로 120ms 전체 계산을 시작했다.
기존 연속 100Hz·1g 실행과 새 실행의 `circles.csv`, `geometry.vtk`,
`geometry_validation.csv`, `particle_validation.csv`가 SHA-256 기준으로 모두 일치했다.
가진 관련 설정과 모델 식별자를 제외한 저장 설정도 모두 일치했다.

형상 seed: `20261014`.
이번 환경의 `circles.csv` SHA-256:
`90c93b6b90f788e29bc744be6c180b0f26676c760026e9396d8019586827dd95`.

시작 구간의 결과 갱신 및 수치 유한성을 확인했다.
전체 계산 완료, 장시간 안정성, 최종 충전율은 아직 검증하지 않았다.
실행 결과·로그·이전 보관 데이터는 Git에 포함하지 않는다.

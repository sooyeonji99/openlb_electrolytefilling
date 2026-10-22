# 단시간 검증 결과

전체1.06초 계산은 시작하지 않았다.

|검사|결과|범위|
|---|---|---|
|전체 파형 수식|PASS|30ms 대기, 10→100→10Hz 연속 위상, 가속도 미분, 총55주기, 종료30ms 대기|
|구조·입자|PASS|높이3.67mm, 좌우하단 벽·상부 개방, 기공 연결성 및 입자 간격|
|smoke_initial_rest_final|PASS|200step, 4 MPI프로세스, pressure2500Pa·유한성·VTK 검사|
|smoke_active_up_final|PASS|200step, 4 MPI프로세스, pressure2500Pa·유한성·VTK 검사|
|smoke_all_transitions_final|PASS|640step, 4 MPI프로세스, pressure2500Pa·유한성·VTK 검사|

smoke_initial_rest_final은 실제 기본 조건의 초기200step이다. smoke_active_up_final은 pre-rest만0으로 바꾸어 실제 진폭7.34mm와 실제0.5초 상승률의 스윕 시작200step을 검사했다.

smoke_all_transitions_final은 전환을 짧은 계산으로 검사하기 위해 대기5µs·상승15µs·하강15µs·대기5µs, 진폭비0.0002로 바꾼 시험640step이다. 본 실행 설정은 바꾸지 않았다.

초기 시험에서 pre-rest=0일 때 step0 가속도 필드가 중력만 담긴 문제를 발견했다. applyAcceleration(0)을 추가한 후 최종 실행파일로 위 세 검사를 재실행했다. 기존 smoke_active_up은 수정 전 기록이다.

이 검증은 전체1.06초의 장시간 안정성·시간/격자 수렴·상부 액체 도달을 입증하지 않는다. 상부 배출은 개방 경계로 허용되며, 짧은 시험에서 액체가 상부에 도달한 것은 아니다.

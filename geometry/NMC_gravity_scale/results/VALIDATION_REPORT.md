# NMC_gravity_scale 단시간 검증

|검사|판정|범위|
|---|---|---|
|smoke_0g|PASS|4프로세스MPI,200step,가속도/압력/유한성/VTK 검사|
|smoke_1g|PASS|4프로세스MPI,200step,가속도/압력/유한성/VTK 검사|
|smoke_2g|PASS|4프로세스MPI,200step,가속도/압력/유한성/VTK 검사|
|처방 파형|PASS|0·0.5·1·2g,100Hz12주기,분기/위상/가속도 크기|

기본조건 120ms=1,920,000step이다. 짧은MPI검증은200step=12.5µs의초기구간이며, 120ms의유동안정성이나상부액체도달을검증한것은아니다.

M=0은이전모델의무가진초기구간과대조했다. 전체격자·기공 연결성·입자간격도geometry_validation과particle_validation에서검사했다. 장시간계산은자동시작하지않았다.

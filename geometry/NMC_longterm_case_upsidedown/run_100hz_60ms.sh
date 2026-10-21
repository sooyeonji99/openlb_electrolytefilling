#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
exec ./run_mpi.sh --ranks "${1:-4}" \
  --outer-boundary open \
  --dx-um 5 \
  --steps 960000 \
  --electrolyte-pressure-pa 0 \
  --gravity-y -9.81 \
  --excitation-frequency-hz 100 --excitation-height-factor 2 \
  --excitation-start-s 0 --excitation-duration-s 0.06 \
  --excitation-on-s 0 --excitation-off-s 0 \
  --tau 0.5075 --phase-tau 0.5075 \
  --output-every 10000 --diagnostic-every 1 \
  --output "${2:-results/upsidedown_100hz_60ms_$(date +%Y%m%d_%H%M%S)}"

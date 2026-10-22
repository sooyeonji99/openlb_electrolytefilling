#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
exec ./run_mpi.sh --ranks "${1:-4}" \
  --outer-boundary top-open --dx-um 5 \
  --duration-s 1.06 \
  --electrolyte-pressure-pa 2500 --pressure-ramp-s 0 \
  --gravity-y -9.81 \
  --pre-rest-s 0.03 --sweep-up-s 0.5 --sweep-down-s 0.5 --post-rest-s 0.03 \
  --sweep-low-hz 10 --sweep-high-hz 100 --excitation-height-factor 2 \
  --tau 0.5075 --phase-tau 0.5075 \
  --output-every 10000 --diagnostic-every 1 \
  --output "${2:-results/forcesweep_$(date +%Y%m%d_%H%M%S)}"

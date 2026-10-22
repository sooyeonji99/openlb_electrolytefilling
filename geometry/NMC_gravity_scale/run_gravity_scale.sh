#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
scale="${1:-1}"
output="${2:-results/gravity_scale_${scale}g_$(date +%Y%m%d_%H%M%S)}"
exec ./run_mpi.sh --ranks "${3:-4}" \
  --outer-boundary top-open --dx-um 5 \
  --duration-s 0.12 \
  --electrolyte-pressure-pa 2500 --pressure-ramp-s 0 \
  --gravity-y -9.81 \
  --excitation-frequency-hz 100 --gravity-scale "$scale" \
  --tau 0.5075 --phase-tau 0.5075 \
  --output-every 10000 --diagnostic-every 1 \
  --output "$output"

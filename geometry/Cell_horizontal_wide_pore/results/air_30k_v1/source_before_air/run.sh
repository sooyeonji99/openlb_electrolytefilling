#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
steps="${1:-2000}"
if [[ ! "$steps" =~ ^[1-9][0-9]*$ ]]; then
  printf 'Step count must be a positive integer\n' >&2
  exit 1
fi
run_dir="results/wide_pore_theta30_p2500_${steps}_$(date +%Y%m%d_%H%M%S)_$$"
printf 'Output: %s\n' "$run_dir"
./main --contact-angle-deg 30 --dx-um 5 --tau 0.53 --phase-tau 0.53 \
  --nu-liquid 1e-6 --nu-gas 1e-7 \
  --outer-boundary open --ns-collision bgk --source-mode legacy \
  --reservoir-pressure-pa 2500 --pressure-ramp-s 0.00025 \
  --steps "$steps" --diagnostic-every 1 --output-every 100 \
  --output "$run_dir" > "$run_dir.log" 2>&1
python3 verify_vtk.py "$run_dir" > "$run_dir.vtk_check.log" 2>&1
printf 'Completed: %s\n' "$run_dir"

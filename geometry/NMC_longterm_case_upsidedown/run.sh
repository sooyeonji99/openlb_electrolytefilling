#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
if [[ "${1:-}" == --help || "${1:-}" == -h ]]; then
  exec ./main --help
fi
args=("$@")
run_dir="results/run_$(date +%Y%m%d_%H%M%S)_$$"
for ((i=0; i<${#args[@]}; i+=2)); do
  if [[ "${args[i]}" != --* || $((i+1)) -ge ${#args[@]} ]]; then
    printf 'Options require --name value pairs; use --help.\n' >&2
    exit 2
  fi
  if [[ "${args[i]}" == --output ]]; then run_dir="${args[i+1]}"; fi
done
if [[ -e "$run_dir" || -e "$run_dir.log" ]]; then
  printf 'Output already exists: %s. Choose a new --output directory.\n' "$run_dir" >&2
  exit 2
fi
mkdir -p -- "$(dirname -- "$run_dir")"
printf 'Output: %s\nLog: %s.log\n' "$run_dir" "$run_dir"
status=0
./main --output "$run_dir" "${args[@]}" > "$run_dir.log" 2>&1 || status=$?
if [[ -f "$run_dir/multiphase_diagnostics.csv" ]]; then
  python3 analyze_runs.py "$run_dir" > "$run_dir.analysis.log" 2>&1 || printf 'Analysis needs inspection: %s.analysis.log\n' "$run_dir" >&2
fi
if [[ "$status" -eq 0 ]]; then
  python3 verify_vtk.py "$run_dir" > "$run_dir.vtk_check.log" 2>&1
  python3 verify_upsidedown_layout.py "$run_dir" > "$run_dir.orientation_check.log" 2>&1
  printf 'Completed: %s\n' "$run_dir"
else
  printf 'Simulation stopped (exit %s). Inspect %s.log and analysis_summary.json\n' "$status" "$run_dir" >&2
  exit "$status"
fi

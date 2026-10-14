#!/usr/bin/env python3
import csv
from pathlib import Path

root = Path("results/geometry_validation")
geometry = next(csv.DictReader((root / "geometry_validation.csv").open(newline="")))
result = next(csv.DictReader((root / "validation_result.csv").open(newline="")))
costs = list(csv.DictReader((root / "time_cost.csv").open(newline="")))
multiphase_path = Path("results/multiphase_smoke/multiphase_validation.csv")
multiphase = next(csv.DictReader(multiphase_path.open(newline="")))
dx2 = next(row for row in costs if abs(float(row["dx_um"])-2.0) < 1e-12)
print(f"status={result['status']}")
print(f"outer_area_error={float(geometry['outer_area_relative_error']):.6%}")
print(f"electrolyte_area_error={float(geometry['electrolyte_area_relative_error']):.6%}")
print(f"annulus_porosity={float(geometry['annulus_porosity']):.6f}")
print(f"connected_pore_fraction={float(geometry['connected_pore_fraction']):.6f}")
print(f"dx2_cells={float(dx2['quarter_sector_cells']):.0f}")
print(f"dx2_memory_GiB={float(dx2['estimated_memory_GiB']):.3f}")
print(f"dx2_steps_10s={float(dx2['steps_requested']):.0f}")
print(f"multiphase_status={multiphase['status']}")
print(f"multiphase_phi_range=[{multiphase['min_phi']},{multiphase['max_phi']}]")
print(f"multiphase_max_speed_lattice={multiphase['max_speed_lattice']}")
raise SystemExit(0 if result["status"] == "PASS" and multiphase["status"] == "PASS" else 1)

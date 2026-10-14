#!/usr/bin/env python3
import csv
import math
from pathlib import Path

rows = list(csv.DictReader(Path("results/integrated/small_pore_diagnostics.csv").open(newline="")))
final = rows[-1]
front = float(final["front_m"])
exit_position = float(final["second_throat_exit_m"])
saturation = float(final["liquid_saturation"])
speed = float(final["max_speed_lattice"])
finite = all(math.isfinite(value) for value in (front, saturation, speed))
passed = finite and front >= exit_position and 0 <= saturation <= 1 and speed < 0.01
result = {
    "status": "PASS" if passed else "FAIL",
    "physical_time_s": final["time_s"],
    "front_m": front,
    "second_throat_exit_m": exit_position,
    "front_beyond_second_throat_m": front - exit_position,
    "liquid_saturation": saturation,
    "max_speed_lattice": speed,
}
with Path("small_pore_validation.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=result.keys())
    writer.writeheader(); writer.writerow(result)
print(", ".join(f"{key}={value}" for key, value in result.items()))
raise SystemExit(0 if passed else 1)

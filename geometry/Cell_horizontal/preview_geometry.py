"""Create a vector preview from the shortened geometry validation output."""
import csv
from pathlib import Path
root=Path(__file__).parent
folder=root/'results/short_geometry_dx02'
rows=list(csv.DictReader((folder/'circles.csv').open()))
config=next(csv.DictReader((folder/'geometry_validation.csv').open()))
w=float(config['width_m'])*1e6;h=float(config['height_m'])*1e6
with (root/'geometry_preview.svg').open('w') as f:
    f.write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="-150 -180 {w+300} {h+360}"><rect width="{w}" height="{h}" fill="#e7eef2"/>')
    for c in rows:
        if c['removed_for_source']=='0':
            f.write(f'<circle cx="{float(c["x_m"])*1e6:.2f}" cy="{h-float(c["y_m"])*1e6:.2f}" r="{float(c["radius_m"])*1e6:.2f}" fill="#454d58"/>')
    f.write(f'<rect x="{(w-1000)/2}" y="{h-500}" width="1000" height="500" fill="#36a9cf"/><rect width="{w}" height="{h}" fill="none" stroke="#111" stroke-width="20" stroke-dasharray="50 35"/><line x1="0" y1="{h}" x2="{w}" y2="{h}" stroke="#111" stroke-width="20"/><text x="{w/2}" y="-65" font-size="90" text-anchor="middle">4.3 × {h/1000:.3f} mm; source 1 × 0.5 mm</text></svg>')

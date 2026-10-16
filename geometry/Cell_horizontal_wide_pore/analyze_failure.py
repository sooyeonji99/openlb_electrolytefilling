"""Locate phase/speed extrema in the final 1000 steps of an isolated run."""
import array,base64,csv,json,math,re,sys,xml.etree.ElementTree as ET
from pathlib import Path
p=Path(sys.argv[1]);config=dict(csv.reader((p/'run_config.csv').open()))
dx=float(config['dx_m']);dt=float(config['dt_s']);nx_phys=round(float(config['width_m'])/dx);ny_phys=round(float(config['height_m'])/dx)
rows=list(csv.DictReader((p/'multiphase_diagnostics.csv').open()));last=int(rows[-1]['step'])
circles=[{k:float(v) for k,v in r.items()} for r in csv.DictReader((p/'circles.csv').open()) if float(r['removed_for_source'])==0]
result={'last_diagnostic':rows[-1],'snapshots':[]}
for limit in (.1,.3,1):result[f'first_speed_above_{limit}']=next((r for r in rows if float(r['max_speed_lattice'])>limit),None)
for path in sorted((p/'vtkData/data').glob('cell_horizontal_*.vti')):
 step=int(re.search(r'iT(\d+)',path.name)[1])
 if step<last-1000:continue
 root=ET.parse(path);ext=list(map(int,root.find('.//Piece').get('Extent').split()));nx=ext[1]-ext[0]+1;f={}
 for n in root.findall('.//DataArray'):
  a=array.array('f');a.frombytes(base64.b64decode(''.join(n.text.split())[8:]));f[n.get('Name')]=a
 mat=f['material'];phi=f['phi_fluid_only'];vel=f['velocity_phys_m_s'];pressure=f['pressure_phys_Pa'];best=(-1,None);lo=(float('inf'),None);hi=(-float('inf'),None);outlet=0
 for i,m in enumerate(mat):
  if m not in (1,3,4,5,7):continue
  x=i%nx+ext[0];y=i//nx+ext[2]
  if not(0<=x<=nx_phys and 0<=y<=ny_phys):continue
  u=math.hypot(vel[2*i],vel[2*i+1])*dt/dx
  if u>best[0]:best=(u,i)
  if phi[i]<lo[0]:lo=(phi[i],i)
  if phi[i]>hi[0]:hi=(phi[i],i)
  if m in (4,5,7):outlet=max(outlet,u)
 def point(i):
  x=(i%nx+ext[0])*dx;y=(i//nx+ext[2])*dx
  near=sorted((math.hypot(x-c['x_m'],y-c['y_m'])-c['radius_m'],j) for j,c in enumerate(circles))[:2]
  return {'x_um':x*1e6,'y_um':y*1e6,'material':int(mat[i]),'phi':phi[i],'pressure_pa':pressure[i],'nearest_circles':[{'index':j,'surface_distance_um':d*1e6,'center_x_um':circles[j]['x_m']*1e6,'center_y_um':circles[j]['y_m']*1e6,'radius_um':circles[j]['radius_m']*1e6} for d,j in near]}
 r={'step':step,'max_speed_lattice_vtk':best[0],'speed_point':point(best[1]),'min_phi_point':point(lo[1]),'max_phi_point':point(hi[1]),'max_outlet_speed_lattice_vtk':outlet};result['snapshots'].append(r)
 print(step,'U',round(best[0],6),'at',r['speed_point']['x_um'],r['speed_point']['y_um'],'min_phi',lo[0],flush=True)
(p/'failure_analysis.json').write_text(json.dumps(result,indent=2)+'\n')

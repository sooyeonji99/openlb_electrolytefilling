"""Validate physical orientation in the saved initial MPI VTI blocks."""
import array,base64,csv,json,math,sys,xml.etree.ElementTree as ET
from pathlib import Path
p=Path(sys.argv[1]);cfg=dict(csv.reader((p/'run_config.csv').open()));dx=float(cfg['dx_m']);height=float(cfg['height_m']);width=float(cfg['width_m']);cells={};source_phi=[]
for f in (p/'vtkData/data').glob('cell_horizontal_iT0000000iC*.vti'):
 root=ET.parse(f);image=root.find('ImageData');origin=list(map(float,image.attrib['Origin'].split()));extent=list(map(int,image.find('Piece').attrib['Extent'].split()));nx=extent[1]-extent[0]+1
 fields={}
 for node in image.findall('.//DataArray'):
  raw=''.join(node.text.split());v=array.array('f');v.frombytes(base64.b64decode(raw[8:]));fields[node.attrib['Name']]=v
 for i,m in enumerate(fields['material']):
  ix=round((origin[0]+(i%nx+extent[0])*dx)/dx);iy=round((origin[1]+(i//nx+extent[2])*dx)/dx)
  if 0<=ix<=round(width/dx) and 0<=iy<=round(height/dx):cells[ix,iy]=int(m)
  if m==3:
   y=origin[1]+(i//nx+extent[2])*dx
   assert y>=height-float(cfg['source_height_m'])-dx*1e-5
   source_phi.append(fields['phi_fluid_only'][i])
assert source_phi and all(abs(v-1)<1e-6 for v in source_phi)
nx=round(width/dx);ny=round(height/dx)
assert len(cells)==(nx+1)*(ny+1)
assert all(cells[i,ny-1]==2 for i in range(1,nx))
assert all(cells[i,1]==6 for i in range(2,nx-1))
assert all(cells[1,j]==4 and cells[nx-1,j]==5 for j in range(2,ny-1))
assert not any(m==7 for m in cells.values())
result={'status':'PASS','top_source_phi':1,'top_wall':True,'bottom_open':True,'left_right_open':True,'unique_geometry_nodes':len(cells)}
(p/'orientation_validation.json').write_text(json.dumps(result,indent=2));print(result)

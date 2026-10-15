"""Closed-cell phase-field diagnostics; plateau is a candidate, not physical validation."""
import argparse,array,base64,csv,json,math,re
from collections import deque
from pathlib import Path
import xml.etree.ElementTree as ET

parser=argparse.ArgumentParser()
parser.add_argument('run',type=Path)
parser.add_argument('--window-s',type=float,default=1e-4)
parser.add_argument('--absolute-tolerance',type=float,default=1e-6)
parser.add_argument('--relative-tolerance',type=float,default=1e-3)
args=parser.parse_args()
if args.window_s<=0 or args.absolute_tolerance<0 or args.relative_tolerance<0:parser.error('invalid thresholds')
folder=args.run
rows=list(csv.DictReader((folder/'multiphase_diagnostics.csv').open()))
config=dict(csv.reader((folder/'run_config.csv').open()))
if config.get('outer_boundary')=='open':parser.error('Open boundary run: use analyze.py and verify_vtk.py; closed-cell stopping analysis does not apply.')
geometry=next(csv.DictReader((folder/'geometry_validation.csv').open()))
result=next(csv.DictReader((folder/'multiphase_validation.csv').open()))
dx=float(config['dx_m']);pore_cells=int(geometry['pore_cells']);pore_area=pore_cells*dx*dx
stable=result['status']=='PASS' and all(r['finite']=='1' for r in rows)
end=float(rows[-1]['time_s']); windows=[]
for w in range(3):
    lo=end-(3-w)*args.window_s;hi=end-(2-w)*args.window_s
    selected=[r for r in rows if lo-1e-14<=float(r['time_s'])<=hi+1e-14]
    if lo < -1e-14 or len(selected)<3:continue
    values=[float(r['pore_mean_phi']) for r in selected];times=[float(r['time_s']) for r in selected]
    if times[-1]-times[0]<.95*args.window_s:continue
    span=max(values)-min(values);change=values[-1]-values[0]
    tolerance=min(args.absolute_tolerance,args.relative_tolerance*max(abs(values[-1]),1e-12))
    windows.append({'start_s':times[0],'end_s':times[-1],'saturation_span':span,'saturation_change':change,'mean_rate_s_inverse':change/(times[-1]-times[0]),'tolerance':tolerance,'pass':span<=tolerance})
plateau=len(windows)==3 and all(w['pass'] for w in windows)
status='NUMERICAL_FAILURE' if not stable else ('PLATEAU_CANDIDATE' if plateau else 'NOT_YET_DEMONSTRATED')
with (folder/'closed_cell_history.csv').open('w') as f:
    writer=csv.writer(f);writer.writerow(['step','time_s','pore_mean_phi_raw','gas_area_raw_m2','gas_fraction_raw','finite'])
    for r in rows:
        saturation=float(r['pore_mean_phi']);writer.writerow([r['step'],r['time_s'],saturation,(1-saturation)*pore_area,1-saturation,r['finite']])
summary={'geometry_version':config.get('geometry_version','legacy_single_diameter_v1'),'status':status,'numerical_checks_pass':stable,'last_step':int(rows[-1]['step']),'time_s':end,'initial_saturation_raw':float(rows[0]['pore_mean_phi']),'final_saturation_raw':float(rows[-1]['pore_mean_phi']),'gas_fraction_raw':1-float(rows[-1]['pore_mean_phi']),'pore_area_m2':pore_area,'max_sampled_speed_lattice':max(float(r['max_speed_lattice']) for r in rows),'plateau_window_s':args.window_s,'required_windows':3,'absolute_tolerance':args.absolute_tolerance,'relative_tolerance':args.relative_tolerance,'windows':windows,'limitations':['No gas vent; incompressible gas, no dissolution/compression model.','Phase-field area is a 2D proxy, not gas mass.','Plateau candidate requires longer confirmation, conservation and grid/time convergence.']}
# Analyze latest serial snapshot. Never silently join MPI block halos.
files=sorted((folder/'vtkData'/'data').glob('cell_horizontal_*iC*.vti'))
if files:
    latest=max(int(re.search(r'iT(\d+)',p.name)[1]) for p in files)
    selected=[p for p in files if int(re.search(r'iT(\d+)',p.name)[1])==latest]
    if len(selected)!=1:summary['snapshot_note']='MPI topology analysis requires merged blocks; skipped.'
    else:
        root=ET.parse(selected[0]);extent=list(map(int,root.find('.//Piece').attrib['Extent'].split()));nx=extent[1]-extent[0]+1;ny=extent[3]-extent[2]+1
        fields={}
        for node in root.findall('.//DataArray'):
            if node.attrib['Name'] not in ('material','phi_fluid_only'):continue
            values=array.array('f');values.frombytes(base64.b64decode(''.join(node.text.split())[8:]));fields[node.attrib['Name']]=values
        material=fields['material'];phi=fields['phi_fluid_only'];gas=bytearray(nx*ny);liquid=bytearray(nx*ny);seen=bytearray(nx*ny)
        # Exclude the one-cell VTK output halo using physical domain bounds.
        for y in range(ny):
            if not 0<=y+extent[2]<=round(float(config['height_m'])/dx):continue
            for x in range(nx):
                if not 0<=x+extent[0]<=round(float(config['width_m'])/dx):continue
                k=y*nx+x
                if material[k]==1 and math.isfinite(phi[k]):gas[k]=phi[k]<.5;liquid[k]=phi[k]>=.5
                elif material[k]==3:liquid[k]=phi[k]>=.5
        def neighbours(k):
            x=k%nx;y=k//nx
            if x:yield k-1
            if x+1<nx:yield k+1
            if y:yield k-nx
            if y+1<ny:yield k+nx
        components=[]
        for k,occupied in enumerate(gas):
            if not occupied or seen[k]:continue
            queue=deque([k]);seen[k]=1;count=0;sx=sy=0;xmin=ymin=10**9;xmax=ymax=-1
            while queue:
                n=queue.popleft();x=n%nx+extent[0];y=n//nx+extent[2];count+=1;sx+=x;sy+=y;xmin=min(xmin,x);xmax=max(xmax,x);ymin=min(ymin,y);ymax=max(ymax,y)
                for nb in neighbours(n):
                    if gas[nb] and not seen[nb]:seen[nb]=1;queue.append(nb)
            components.append({'nodes':count,'area_m2':count*dx*dx,'centroid_x_m':sx/count*dx,'centroid_y_m':sy/count*dx,'min_x_m':xmin*dx,'max_x_m':xmax*dx,'min_y_m':ymin*dx,'max_y_m':ymax*dx})
        components.sort(key=lambda c:c['nodes'],reverse=True)
        with (folder/'residual_gas_components.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=['nodes','area_m2','centroid_x_m','centroid_y_m','min_x_m','max_x_m','min_y_m','max_y_m']);writer.writeheader();writer.writerows(components)
        queue=deque(k for k,v in enumerate(liquid) if v and material[k]==3);seen=bytearray(nx*ny)
        for k in queue:seen[k]=1
        connected=0
        while queue:
            k=queue.popleft();connected+=material[k]==1
            for nb in neighbours(k):
                if liquid[nb] and not seen[nb]:seen[nb]=1;queue.append(nb)
        summary['snapshot']={'step':latest,'is_final_step':latest==summary['last_step'],'gas_threshold_phi':.5,'gas_nodes':sum(gas),'gas_components_4_connected':len(components),'largest_component_nodes':components[0]['nodes'] if components else 0,'source_connected_liquid_pore_nodes':connected,'source_connected_liquid_pore_fraction':connected/pore_cells,'note':'Threshold components are geometric proxies, not proof of physically trapped bubbles.'}
(folder/'closed_cell_analysis.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
print(json.dumps(summary,indent=2,allow_nan=False))

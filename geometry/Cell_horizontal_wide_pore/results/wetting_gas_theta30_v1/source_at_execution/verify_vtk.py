"""Read OpenLB inline base64 Float32 VTI with the Python standard library."""
import array,base64,csv,math,xml.etree.ElementTree as ET
from pathlib import Path
import sys,re
folder=Path(sys.argv[1])
config=dict(csv.reader((folder/"run_config.csv").open()))
expected_pressure=float(config["pressure_pa"])
checks=[]
for path in sorted((folder/'vtkData'/'data').glob('cell_horizontal_*.vti')):
    step=int(re.search(r'iT(\d+)',path.name)[1]);time=step*float(config['dt_s'])
    ramp=float(config.get('pressure_ramp_s',0))
    expected_pressure=float(config['pressure_pa'])*(min(1,max(0,time/ramp)) if ramp>0 else 1)
    fields={}
    for node in ET.parse(path).findall('.//DataArray'):
        raw=''.join(node.text.split()); values=array.array('f')
        values.frombytes(base64.b64decode(raw[8:]))
        fields[node.attrib['Name']]=values
    material=fields['material'];phi=fields['phi_fluid_only'];pressure=fields['pressure_phys_Pa'];velocity=fields['velocity_phys_m_s']
    source_pressure=[p for m,p in zip(material,pressure) if m==3]
    valid=all(math.isfinite(phi[i]) and math.isfinite(pressure[i]) and math.isfinite(velocity[2*i]) and math.isfinite(velocity[2*i+1]) for i,m in enumerate(material) if m in (1,3,4,5,6,7))
    pressure_ok=config.get('source_mode')=='finite' or all(abs(p-expected_pressure)<.02 for p in source_pressure)
    valid=valid and pressure_ok
    checks.append({'file':path.name,'finite_fluid_fields':int(valid),'source_pressure_min_Pa':min(source_pressure),'source_pressure_max_Pa':max(source_pressure)})
    if config.get('outer_boundary')=='open':
        outlet_pressure=[p for m,p in zip(material,pressure) if m in (4,5,6,7)]
        root=ET.parse(path);extent=list(map(int,root.find('.//Piece').attrib['Extent'].split()));nx=extent[1]-extent[0]+1
        net_flux=liquid_flux=0.;count=0;backflow=0
        for i,m in enumerate(material):
            x=i%nx+extent[0];y=i//nx+extent[2]
            if not (0<=x<=round(float(config['width_m'])/float(config['dx_m'])) and 0<=y<=round(float(config['height_m'])/float(config['dx_m']))):continue
            if m not in (4,5,6,7):continue
            ux,uy=velocity[2*i:2*i+2];un={4:-ux,5:ux,6:-uy,7:uy}[int(m)]
            net_flux+=un*float(config['dx_m']);liquid_flux+=phi[i]*un*float(config['dx_m']);count+=1;backflow+=un < -1e-10
        checks[-1].update(outlet_nodes=count,outlet_pressure_min_Pa=min(outlet_pressure),outlet_pressure_max_Pa=max(outlet_pressure),outlet_net_advective_flux_m2_s=net_flux,outlet_liquid_advective_flux_m2_s=liquid_flux,outlet_backflow_nodes=backflow)
        if not all(math.isfinite(p) and abs(p)<.02 for p in outlet_pressure):checks[-1]['finite_fluid_fields']=0
with (folder/'vtk_validation.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(checks[0]));writer.writeheader();writer.writerows(checks)
print(checks)

if not all(row["finite_fluid_fields"] for row in checks):raise SystemExit("VTK validation failed")

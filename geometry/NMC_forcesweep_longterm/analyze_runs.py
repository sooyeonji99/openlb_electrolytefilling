"""Summarize genuine runs and locate anomaly VTI extrema; no external packages."""
import array,base64,csv,json,math,re,sys,xml.etree.ElementTree as ET
from pathlib import Path

def safe_json(value):
    if isinstance(value,float) and not math.isfinite(value):return str(value)
    if isinstance(value,dict):return {k:safe_json(v) for k,v in value.items()}
    if isinstance(value,list):return [safe_json(v) for v in value]
    return value

def read_vti(path):
    tree=ET.parse(path); im=tree.find('.//ImageData');piece=tree.find('.//Piece')
    extent=list(map(int,piece.attrib['Extent'].split()))
    origin=list(map(float,im.attrib['Origin'].split()));spacing=list(map(float,im.attrib['Spacing'].split()))
    fields={}
    for node in piece.findall('.//DataArray'):
        raw=''.join(node.text.split()); a=array.array('f');a.frombytes(base64.b64decode(raw[8:]));fields[node.attrib['Name']]=a
    def location(i):
        nx=extent[1]-extent[0]+1
        return [origin[0]+(i%nx+extent[0])*spacing[0],origin[1]+(i//nx+extent[2])*spacing[1]]
    return fields,location,extent

def locate(path,kind,circles):
    fields,location,extent=read_vti(path);material=fields['material']
    fluid=[i for i,m in enumerate(material) if m in (1,3,4,5,6,7)]
    if kind=='raw':
        phi=fields['phi_raw_populations'];finite=[i for i in fluid if math.isfinite(phi[i])]
        indices=[]
        if finite:indices += [('min_phi',min(finite,key=lambda i:phi[i])),('max_phi',max(finite,key=lambda i:phi[i]))]
        indices += [('nonfinite_phi',i) for i in fluid if not math.isfinite(phi[i])][:12]
        nonfinite=sum(not math.isfinite(phi[i]) for i in fluid)
    else:
        phi=fields['phi_fluid_only'];v=fields['velocity_phys_m_s'];valid=[i for i in fluid if math.isfinite(v[2*i]) and math.isfinite(v[2*i+1])];bad=[i for i in fluid if not(math.isfinite(v[2*i]) and math.isfinite(v[2*i+1]))];indices=([('max_speed',max(valid,key=lambda i:math.hypot(v[2*i],v[2*i+1])))] if valid else [])+[('nonfinite_velocity',i) for i in bad[:8]];nonfinite=sum(not math.isfinite(phi[i]) for i in fluid)
    points=[]
    nx=extent[1]-extent[0]+1
    for label,i in indices:
        x,y=location(i);nearest=min(circles,key=lambda c:abs(math.hypot(x-c[0],y-c[1])-c[2]))
        neighbors=[]
        for dy in (-1,0,1):
            for dx in (-1,0,1):
                j=i+dy*nx+dx
                if 0<=j<len(material):neighbors.append({'offset':[dx,dy],'material':int(material[j]),'phi':str(phi[j])})
        p={'kind':label,'xy_m':[x,y],'material':int(material[i]),'phi':str(phi[i]),'nearest_circle_surface_distance_m':math.hypot(x-nearest[0],y-nearest[1])-nearest[2],'nearest_circle_xyr_m':nearest,'neighbors':neighbors}
        if kind!='raw':p['velocity_m_s']=[v[2*i],v[2*i+1]];p['speed_m_s']=math.hypot(v[2*i],v[2*i+1])
        if kind!='raw':
            for name in ['rho_lattice','laplace_rho_lattice']:
                if name in fields:p[name]=str(fields[name][i])
            if 'rho_lattice' in fields and 'laplace_rho_lattice' in fields:p['corrected_density_lattice']=str(fields['rho_lattice'][i]+fields['laplace_rho_lattice'][i]/4)
            for name in ['total_acceleration_lattice','phase_force_lattice','external_acceleration_lattice']:
                if name in fields:p[name]=[str(x) for x in fields[name][2*i:2*i+2]]
        points.append(p)
    return {'file':path.name,'nonfinite_fluid_phi_count':nonfinite,'points':points}

def summarize(folder):
    cfg=dict(list(csv.reader((folder/'run_config.csv').open()))[1:])
    count=0;maxu=0;rawmin=math.inf;rawmax=-math.inf;l1=0;residual=0;all_finite=True;first_speed=None;worst_step=0;range_events=0;last=None
    with (folder/'multiphase_diagnostics.csv').open() as f:
        for row in csv.DictReader(f):
            count+=1;last=row;step=int(row['step']);u=float(row['max_speed_lattice']);lo=float(row['raw_min_phi']);hi=float(row['raw_max_phi'])
            maxu=max(maxu,u)
            if lo<rawmin:rawmin=lo;worst_step=step
            rawmax=max(rawmax,hi);l1+=float(row['phase_projection_l1']);residual=max(residual,abs(float(row['phase_projection_mass_residual'])))
            all_finite &= row['finite']=='1'
            if u>.3 and first_speed is None:first_speed=step
            range_events += lo<-.1 or hi>1.1
    if last is None:raise RuntimeError('No complete diagnostic rows')
    log=folder.with_suffix('.log').read_text(errors='replace')
    finalpath=folder/'multiphase_validation.csv';final=list(csv.DictReader(finalpath.open()))[0] if finalpath.exists() else None
    exception=re.search(r'Exception at step (\d+): ([^\n]+)',log)
    requested=int(cfg['requested_steps']);step=int(last['step'])
    status='PASS' if final and final['status']=='PASS' and int(final['steps'])==requested else 'FAIL' if exception or final or 'NMC_longterm_case:' in log else 'INCOMPLETE'
    events=list(csv.DictReader((folder/'anomaly_events.csv').open()))
    forcing_last=None
    if (folder/'forcing_history.csv').exists():
        with (folder/'forcing_history.csv').open() as f:
            for forcing_last in csv.DictReader(f):pass
    summary={'case':folder.name,'status':status,'requested_steps':requested,'last_logged_step':step,'failure_step':int(exception[1]) if exception else None,'reason':exception[2] if exception else None,'last_physical_time_s':float(last['time_s']),'pressure_pa':float(cfg['pressure_pa']),'gravity_y_m_s2':float(cfg['gravity_y_m_s2']),'phase_tau':float(cfg['phase_tau']),'nu_gas_m2_s':float(cfg['nu_gas_m2_s']),'pressure_ramp_s':float(cfg['pressure_ramp_s']),'pressure_duration_s':float(cfg['pressure_duration_s']),'source_mode':cfg['source_mode'],'max_speed_lattice':maxu,'raw_min_phi':rawmin,'raw_max_phi':rawmax,'projection_l1_total':l1,'projection_mass_residual_max':residual,'final_pore_mean_phi':float(last['pore_mean_phi']),'event_rows':len(events),'diagnostic_rows':count,'all_logged_finite':all_finite,'first_speed_over_0_3':first_speed,'raw_range_warning_steps':range_events,'forcing_last':forcing_last,'pressure_release_reached':float(last['time_s'])>=float(cfg['pressure_duration_s']),'excitation_cycles_observed':float(cfg['excitation_frequency_hz'])*max(0.,min(float(last['time_s'])-float(cfg.get('excitation_start_s',0)),float(cfg.get('excitation_duration_s','inf')))),'config':cfg}
    if cfg.get('excitation_mode')=='linear_frequency_sweep':
        pre=float(cfg['sweep_pre_rest_s']);up=float(cfg['sweep_up_s']);down=float(cfg['sweep_down_s'])
        low=float(cfg['sweep_low_hz']);high=float(cfg['sweep_high_hz']);elapsed=max(0.,float(last['time_s'])-pre)
        a=min(elapsed,up);b=min(max(0.,elapsed-up),down)
        summary['excitation_cycles_observed']=low*a+.5*(high-low)/up*a*a+high*b+.5*(low-high)/down*b*b
        summary['excitation_mode']='linear_frequency_sweep'
    outlet_path=folder/'outlet_diagnostics.csv'
    if outlet_path.exists():
        outlet_max=0.;outlet_last=None;first_outlet=None
        with outlet_path.open() as f:
            for row in csv.DictReader(f):
                outlet_last=row;value=float(row['max_outlet_and_inner_speed_lattice']);outlet_max=max(outlet_max,value)
                if first_outlet is None and (value>.1 or float(row['nonfinite_fields'])>0):first_outlet=int(row['step'])
        summary['max_outlet_speed_lattice']=outlet_max;summary['first_outlet_anomaly_step']=first_outlet;summary['last_outlet_diagnostics']=outlet_last
        if status=='FAIL' and not summary['reason']:
            summary['reason']='Nonfinite outlet fields' if outlet_last and float(outlet_last['nonfinite_fields'])>0 else 'Final numerical validity checks failed'
    summary['approx_elapsed_s']=(finalpath if finalpath.exists() else folder/'multiphase_diagnostics.csv').stat().st_mtime-(folder/'run_config.csv').stat().st_mtime
    with (folder/'circles.csv').open() as f:circles=[[float(r['x_m']),float(r['y_m']),float(r['radius_m'])] for r in csv.DictReader(f)]
    selected=[];speeds=[e for e in events if e['reason'] in ('speed_threshold','nonfinite_speed','outlet_velocity_or_nonfinite','ns_state_on_exception')]
    for e in speeds[:1]+speeds[-2:]:selected.append(('regular',int(e['step'])))
    raws=[e for e in events if e['stage']=='pre_projection']
    for e in raws[:1]+raws[-2:]:selected.append(('raw',int(e['step'])))
    selected.append(('raw',worst_step));spatial=[]
    for kind,k in dict.fromkeys(selected):
        prefix='anomaly_raw' if kind=='raw' else 'cell_horizontal';p=folder/f'vtkData/data/{prefix}_iT{k:07d}iC00000.vti'
        if not p.exists():continue
        observation=locate(p,kind,circles)
        fields,xy,ext=read_vti(p)
        if kind=='raw':
            bad=[i for i,(v,m) in enumerate(zip(fields['phi_raw_populations'],fields['material'])) if m in (1,3,4,5,6,7) and not math.isfinite(v)]
            if bad:observation['nonfinite_bbox_m']=[[min(xy(i)[j] for i in bad),max(xy(i)[j] for i in bad)] for j in (0,1)]
        spatial.append(observation)
    summary['spatial']=spatial
    (folder/'analysis_summary.json').write_text(json.dumps(safe_json(summary),indent=2,allow_nan=False)+'\n')
    return summary

if __name__=='__main__':
    root=Path(sys.argv[1]);out=[]
    if (root/'run_config.csv').exists():
        result=summarize(root);print(result['case'],result['status'],result['last_logged_step'],result['reason']);raise SystemExit(0)
    for p in sorted(root.iterdir()):
        if p.is_dir() and (p/'run_config.csv').exists():
            result=summarize(p);out.append(result);print(result['case'],result['status'],result['last_logged_step'],result['reason'])
    (root/'comparison.json').write_text(json.dumps(safe_json(out),indent=2,allow_nan=False)+'\n')

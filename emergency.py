import numpy as np, full_ship as fs, warnings; warnings.filterwarnings('ignore')
from contour_model import drag; from scale_model import sail_mass
c=fs.c; ly=fs.ly; yr=fs.yr; sp=fs.sp_eff(**fs.ENG['central'])
r=fs.optimise(100,sp=sp); vc=r['vc']; M=r['M_dry']; rem=4.24-r['s1']
v=np.geomspace(vc,1e5,600)
for budget in (0.5,1,2,5):
    rows=[]
    for R in np.geomspace(1e6,2e8,80):
        for I in np.geomspace(1e5,3e7,80):
            F,_,_=drag(v,R,I,0.06); mw,ms,mt=sail_mass(R,I,Fmax=F.max()); Ms=mw+ms+mt
            if Ms/1e9>budget: continue
            inv=(M+Ms)/F; s=np.trapezoid((inv*v)[::-1],v[::-1]); t=np.trapezoid(inv[::-1],v[::-1])
            rows.append((s,t,R,I,mw/Ms,ms/Ms,mt/Ms))
    rows=np.array(rows); j=np.argmin(rows[:,0]); smin=rows[j,0]; near=rows[rows[:,0]<=1.02*smin]
    print(f"{budget} Mt: {smin/ly:.2f} ly; time {near[:,1].min()/yr:.0f}-{near[:,1].max()/yr:.0f} yr; R {rows[j,2]/1e3:.0f} km I {rows[j,3]/1e6:.2f} MA; wire/hoop/tether {rows[j,4]:.2f}/{rows[j,5]:.2f}/{rows[j,6]:.2f}; {'stops short' if smin/ly<rem else 'overshoots'} (margin {rem-smin/ly:+.2f} ly)",flush=True)

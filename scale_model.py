"""sail_mass(): wire, hoop structure and tether mass (used by the paper).
NOTE: drag() below is the older model, kept for reference only; the paper uses contour_model.drag().
Drag law: Fujita (2004) hybrid-PIC fit, F = Cd (rho v^2/2) pi L^2, Cd = 2.9 fluid limit,
with Fujita kinetic correction using r_iL ~ d_i.
Standoff (Fujita convention): L = (mu0 I^2 R^4 / (16 rho v^2))^(1/6).
If L < R (dipole invalid) -> 'doughnut' tube around the wire: d = mu0 I/(2 pi B_mp), B_mp = sqrt(mu0 rho v^2),
drag area = 2 pi R * 2d (capped at pi (R+d)^2), same Cd and kinetic correction with L->d.
Mass: REBCO wire + hoop-stress structure (Zylon) + tethers to ship.
"""
import numpy as np
mu0=4e-7*np.pi; mp=1.6726e-27; c=2.998e8; yr=3.156e7; ly=9.461e15
M_ship=1e9            # 1 Mt dry
v0=0.05*c
ve=8.0e6              # fusion exhaust

def cd_fujita(x):
    x=np.asarray(x,float)
    return np.where(x<1, 2.9*np.exp(-0.35*x**2), 2.4/np.maximum(x,1e-30)*np.exp(-0.17/np.maximum(x,1e-30)**2))

def drag(v,R,I,n_cm3,coef='fujita'):
    n=n_cm3*1e6; rho=n*mp*1.6           # +60% mass from ionised He
    di=2.28e5/np.sqrt(n_cm3)            # proton inertial length, m
    L=(mu0*I**2*R**4/(16*rho*v**2))**(1/6)
    Bmp=np.sqrt(mu0*rho*v**2)
    d=mu0*I/(2*np.pi*Bmp)
    q=0.5*rho*v**2
    dip = L>=R
    area=np.where(dip, np.pi*L**2, np.minimum(4*np.pi*R*d, np.pi*(R+d)**2))
    Lc=np.where(dip, L, d)
    cd=cd_fujita(di/Lc)
    if coef=='freeland':  cd=cd*0.345/0.575       # Freeland, no /2
    if coef=='freeland2': cd=cd*0.217/0.575       # Freeland with /2
    return cd*q*area, dip, di/Lc

def sail_mass(R,I,J=3e9,rho_w=8500,sig=5.8e9/2,rho_s=1560,Fmax=0):
    Aw=I/J; a=np.sqrt(Aw/np.pi)
    m_w=2*np.pi*R*Aw*rho_w
    T=mu0*I**2/(4*np.pi)*(np.log(8*R/a)-1)
    m_s=2*np.pi*R*(T/sig)*rho_s
    m_t=Fmax*np.sqrt(2)*R/(sig/rho_s)
    return m_w, m_s, m_t

def brake(R,I,n,vf,coef='fujita',J=3e9):
    v=np.geomspace(v0,vf,400)
    F,dip,x=drag(v,R,I,n,coef)
    mw,ms,mt=sail_mass(R,I,J=J,Fmax=F.max())
    m=M_ship+mw+ms+mt
    fuel_end=m*(np.exp(vf/ve)-1)      # fusion does the last vf
    mtot=m+fuel_end
    inv=mtot/F
    t=np.trapezoid(inv[::-1],v[::-1])
    s=np.trapezoid((inv*v)[::-1],v[::-1])
    return dict(R=R,I=I,t=t/yr,s=s/ly,m_sail=mw+ms+mt,mw=mw,ms=ms,mt=mt,fuel=fuel_end,
                arr=mtot,dep=mtot*np.exp(v0/ve),dip_start=bool(dip[0]),x0=x[0],xend=x[-1],F0=F.max())

def best(n,tmax,coef='fujita',J=3e9):
    bestr=None
    for R in np.geomspace(1e5,2e8,60):
        for I in np.geomspace(1e4,1e9,80):
            for vf in (2e5,5e5,1e6,2e6):
                r=brake(R,I,n,vf,coef,J)
                if r['t']<=tmax and r['s']<4.0 and (bestr is None or r['dep']<bestr['dep']):
                    r['vf']=vf; bestr=r
    return bestr

def show(tag,r):
    if r is None: print(f"{tag:38s} no solution"); return
    print(f"{tag:38s} R={r['R']/1e3:9.0f} km  I={r['I']/1e6:7.2f} MA  t={r['t']:5.1f} yr  dist={r['s']:.2f} ly  "
          f"sail={r['m_sail']/1e9:6.3f} Mt (wire {r['mw']/1e9:.3f}, hoop {r['ms']/1e9:.3f}, teth {r['mt']/1e9:.3f})  "
          f"vf={r['vf']/1e3:.0f} km/s  arrive={r['arr']/1e9:5.2f} Mt  DEPART={r['dep']/1e9:5.1f} Mt  "
          f"{'dipole' if r['dip_start'] else 'doughnut'} di/L={r['x0']:.3f}")

if __name__=='__main__':
    print("All-fusion reference (accelerate+brake): depart =", round(M_ship*np.exp(2*v0/ve)/1e9,1),"Mt")
    print("Fusion accel + perfect free brake:", round(M_ship*np.exp(v0/ve)/1e9,2),"Mt\n")
    for tmax in (20,30,40,60):
        print(f"--- braking time cap {tmax} yr ---")
        for n in (0.03,0.06,0.12):
            show(f"n={n} Fujita Cd2.9 J=3e9",best(n,tmax))
        show("n=0.06 Freeland(no /2)",best(0.06,tmax,'freeland'))
        show("n=0.06 Freeland(with /2)",best(0.06,tmax,'freeland2'))
        show("n=0.06 Fujita, poor wire J=1e9",best(0.06,tmax,'fujita',1e9))
        show("n=0.06 Fujita, good wire J=5e9",best(0.06,tmax,'fujita',5e9))
        print()

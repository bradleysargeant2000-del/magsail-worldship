"""Independent re-check: time-stepped simulation, written separately from full_ship.py."""
import numpy as np, warnings; warnings.filterwarnings('ignore')
import full_ship as fs
c=2.998e8; yr=3.156e7; ly=9.461e15; mu0=4e-7*np.pi; mp=1.6726e-27
sp=fs.sp_eff(**fs.ENG['central'])
print("A. HABITAT recomputed by hand")
R,r=224.,20.; A=4*np.pi**2*R*r; shield=A*6.5e3
lam=shield/(2*np.pi*R); v=np.sqrt(9.81*R); T=lam*v**2; ring=1.5*T/300e6*7850*2*np.pi*R
print(f"  hull area {A:,.0f} m2, shield {shield/1e9:.3f} Mt, ring tension {T:.2e} N, ring {ring/1e9:.4f} Mt (model: {fs.H['Spin ring + spokes (steel, x1.5 margin)']/1e6:.4f} Mt)")
print(f"  model habitat {fs.M_hab/1e9:.4f} Mt")

print("\nB. MISSION by explicit time-stepping (100-yr central optimum)")
best=fs.optimise(100,sp=sp); a0,vc=best['a0'],best['vc']
print(f"  optimum from model: a0={a0:.4f}, vc={vc/c:.4f}c, depart {best['M_dep']/1e9:.2f} Mt, T {best['T']:.2f} yr")
k=fs.k_tank; ve=fs.ve
M=best['M_dep']; F=a0*M; Meng=0.5*F*ve/sp; Mdry=fs.M_hab+Meng
fuel=(M-Mdry)/(1+k)
dt=3600*24.0; t=0; x=0; vel=0; phase=1; m=M
fuel_rate=F/ve
while True:
    if phase==1:
        dm=fuel_rate*dt; m-=dm*(1+k); vel+=F/m*dt
        if vel>=vc:
            phase=2; t1=t; m_cruise=m
            mm=m; vv=vel; xb=0
            while vv>0: mm-=fuel_rate*dt*(1+k); vv-=F/mm*dt; xb+=vv*dt
            x_brake=fs.D_target-xb
    elif phase==2:
        # brake so we stop at target: start when remaining distance = braking distance (computed analytically below via look-ahead)
        if x>=x_brake: phase=3; t2s=t
    elif phase==3:
        dm=fuel_rate*dt; m-=dm*(1+k); vel-=F/m*dt
        if vel<=-fs.dv_dock*0: 
            break
    x+=vel*dt; t+=dt
print(f"  time-stepped: accel {t1/yr:.2f} yr, cruise mass {m_cruise/1e9:.2f} Mt (model {best['M_cruise']/1e9:.2f}), total {t/yr:.2f} yr (model {best['T']:.2f}), arrive {m/1e9:.3f} Mt (model dry {Mdry/1e9:.3f} + docking fuel), distance {x/ly:.3f} ly")

print("\nC. LOOP HOOP TENSION re-derived from inductance: L=mu0 R(ln(8R/a)-2) -> T = I^2/(4pi) dL/dR")
I=1e6; Rl=2e7; a=np.sqrt(I/3e9/np.pi); h=1.0
Lf=lambda RR: mu0*RR*(np.log(8*RR/a)-2)
T_num=I**2/(4*np.pi)*(Lf(Rl+h)-Lf(Rl-h))/(2*h); T_formula=mu0*I**2/(4*np.pi)*(np.log(8*Rl/a)-1)
print(f"  numeric {T_num:.4e} N vs formula {T_formula:.4e} N")

print("\nD. CONTOUR MODEL continuity across the dipole/doughnut switch")
from contour_model import drag
n=0.06; rho=n*1e6*mp*1.6; vv=0.048*c
for I in (1e6,3e6,1e7):
    Rt=(mu0*I**2/(16*rho*vv**2))**0.5   # radius where L=R
    Fa,_,_=drag(np.array([vv]),Rt*1.001,I,n); Fb,_,_=drag(np.array([vv]),Rt*0.999,I,n)
    print(f"  I={I/1e6:.0f} MA: switch at R={Rt/1e3:,.0f} km, force just inside {Fa[0]:.3e} N vs outside {Fb[0]:.3e} N (jump x{Fa[0]/Fb[0]:.2f})")

print("\nE. LITERATURE REPRODUCTION")
for nn in (0.05,0.1):
    Fi,_,x_=drag(np.array([0.05*c]),2.6e5,7.8e6,nn)
    print(f"  Freeland/Icarus 260 km, 7.8 MA, 0.05c, n={nn}: our force {Fi[0]:.2e} N (Freeland ~1e5 N)")
# Gros 2017: 1,500 t sail, R=1,600 km, from 0.1c; stopping distance 0.37 ly (his model)
def drag_gros(v,R,I,n):
    """Gros (2017) single-particle drag (protons only after /1.6 below)."""
    Ic=1.55e6; lg=np.log(np.maximum(I/((v/c)*Ic),1.0))
    return 2*mp*1.6*n*1e6*v**2*0.081*np.pi*R**2*lg**3
Ms=1.5e6; Rg=1.6e6; vg0=0.1*c; ng=0.1
Ig=0.1*4.2e6
v_=np.geomspace(vg0,1e4,3000); Fg=drag_gros(v_,Rg,Ig,ng)/1.6   # Gros uses protons only
s=np.trapezoid((Ms*v_/Fg)[::-1],v_[::-1])
print(f"  Gros example (his formula, n=0.1, I={Ig/1e6:.2f} MA fixed): stopping distance {s/ly:.2f} ly (Gros: 0.37 ly with optimal varying current)")

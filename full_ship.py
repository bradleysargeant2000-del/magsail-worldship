"""Stardust World-Ship: full end-to-end model (one consistent set of assumptions).
Every input is listed in INPUTS with its source/reasoning. Run: python3 full_ship.py
"""
import numpy as np
c=2.998e8; yr=3.156e7; ly=9.461e15; sigma=5.67e-8; g=9.81; MeV=1.602e-13; u=1.6605e-27

# ---------------- INPUTS ----------------
crew=1000
R_spin=224.0          # m, ring radius: 1 g at 2 rpm
r_tube=20.0           # m, torus tube radius (3 decks inside)
shield_areal=6.5      # t/m^2 water: Globus et al. (NSS) recommend 6-7 t/m2 for deep space (20 mSv/yr)
steel_sig=300e6; steel_rho=7850
ve=8.0e6              # m/s fusion exhaust (Charter assumption; tested 5,000-10,000 km/s)
k_tank=0.02           # tank mass per unit fuel mass (dropped as emptied)
sp_engine=1.0e7       # W/kg engine specific power incl. radiators (Daedalus achieved 2.4e7 with D-He3; D-D has more neutrons -> lower)
D_target=4.24*ly
v_target=30e3         # m/s: match Alpha Cen system's motion (~22 km/s radial, ~32 km/s total)
D_frac_water=3.46e-5
dv_dock=50e3          # m/s: match Alpha Cen motion (~32 km/s) + orbit insertion margin
# engine = core (magnets, driver, nozzle) + radiators for heat absorbed by the ship
eta_c=0.30            # fraction of fusion energy that becomes jet energy (half the 62% charged share)
ENG={'optimistic':dict(sp_core=2.4e7,f_abs=0.01,rad_kgm2=0.5,T_rad=1800),
     'central':   dict(sp_core=2.4e7,f_abs=0.03,rad_kgm2=1.0,T_rad=1500),
     'pessimistic':dict(sp_core=1.0e7,f_abs=0.10,rad_kgm2=3.0,T_rad=1200)}
def sp_eff(sp_core,f_abs,rad_kgm2,T_rad):
    q=2*0.9*sigma*T_rad**4          # W per m2 of two-sided panel
    return 1/(1/sp_core + f_abs*rad_kgm2/(eta_c*q))  # kg D per kg water (156 ppm D/H)

# ---------------- HABITAT ----------------
w=np.sqrt(g/R_spin); rpm=w*60/(2*np.pi); vrim=w*R_spin
A_hull=4*np.pi**2*R_spin*r_tube; V=2*np.pi**2*R_spin*r_tube**2
floor=2*np.pi*R_spin*(34+38+34); farm=60*crew
H={}
H['Radiation shield (water 6.5 t/m2)']=A_hull*shield_areal
H['Pressure hull (steel, 1 atm)']=A_hull*(1e5*r_tube/steel_sig)*steel_rho/1e3
m_line=H['Radiation shield (water 6.5 t/m2)']*1e3/(2*np.pi*R_spin)
H['Spin ring + spokes (steel, x1.5 margin)']=1.5*m_line*vrim**2/steel_sig*steel_rho*2*np.pi*R_spin/1e3
H['Decks, homes, buildings (0.2 t/m2)']=floor*0.2
H['Farms: beds + water (60 m2/person, 0.3 t/m2)']=farm*0.3
H['Air (1 atm)']=V*1.2/1e3
H['Life support, recycling, ship power + radiators']=20000
H['Workshops, spares, raw stock for a century']=100000
H['People + belongings (5 t each)']=crew*5.07
H['Forward dust/gas shield']=10000
M_hab=sum(H.values())*1e3   # kg

# ---------------- FUSION FUEL PHYSICS ----------------
# catalysed D-D: 6 D -> 2 He4 + 2 p + 2 n + 43.2 MeV ; neutrons carry 2.45+14.1 = 16.55 MeV
eps=43.2*MeV/(6*2.0141*u*c**2)          # fraction of fuel rest mass released
E_per_kg=eps*c**2
f_neutron=16.55/43.2
ve_max=c*np.sqrt(2*eps)                 # if ALL energy went into the jet
eta_needed=(ve/ve_max)**2               # fraction of fusion energy that must end up as jet energy

# ---------------- MISSION ----------------
def mission(a0, v_cruise, sp=sp_engine, ve=ve, k=k_tank, M_hab=M_hab):
    """Constant-thrust engine sized for a0 at departure. Same engine brakes. Returns dict."""
    ve_eff=ve/(1+k)                      # tanks are dropped mass that gives no thrust
    dv1=v_cruise; dv2=v_cruise-v_target if v_cruise>v_target else 0
    Rm=np.exp((v_cruise+v_cruise+dv_dock)/ve_eff)   # speed up, slow down, + docking/orbit margin
    den=1-a0*ve*Rm/(2*sp)
    if den<=0: return None
    M_dep=M_hab*Rm/den
    F=a0*M_dep; P_jet=0.5*F*ve; M_eng=P_jet/sp
    M_dry=M_hab+M_eng
    M_cruise=M_dry*np.exp(v_cruise/ve_eff)
    mdot=F/ve*(1+k)                      # total mass loss rate (fuel + dropped tanks)
    t1=(M_dep-M_cruise)/mdot; t2=(M_cruise-M_dry)/mdot
    # distances (exact constant-thrust rocket integrals)
    def burn_dist(M0,M1,mdot,v_start,sign):
        t=np.linspace(0,(M0-M1)/mdot,4000); m=M0-mdot*t
        v=v_start+sign*ve_eff*np.log(M0/m)
        return np.trapezoid(v,t)
    s1=burn_dist(M_dep,M_cruise,mdot,0,+1)
    s2=burn_dist(M_cruise,M_dry,mdot,v_cruise,-1)
    s_cr=D_target-s1-s2
    if s_cr<0: return None
    t_cr=s_cr/v_cruise
    fuel=(M_dep-M_dry)/(1+k); tanks=fuel*k
    P_fusion=P_jet/eta_c
    return dict(a0=a0,vc=v_cruise,M_dep=M_dep,M_dry=M_dry,M_eng=M_eng,M_cruise=M_cruise,fuel=fuel,tanks=tanks,
                F=F,P_jet=P_jet,P_fus=P_fusion,t1=t1/yr,t2=t2/yr,t_cr=t_cr/yr,T=(t1+t2+t_cr)/yr,
                s1=s1/ly,s2=s2/ly,a_brake_end=F/M_dry,water=fuel/D_frac_water)

def optimise(Tmax,**kw):
    best=None
    for a0 in np.geomspace(0.003,0.3,80):
        for vc in np.linspace(0.01,0.08,141)*c:
            r=mission(a0,vc,**kw)
            if r and r['T']<=Tmax and (best is None or r['M_dep']<best['M_dep']): best=r
    return best

if __name__=='__main__':
    print("=== HABITAT ===")
    print(f"spin {rpm:.2f} rpm (1 g at {R_spin:.0f} m), rim speed {vrim:.0f} m/s; volume {V:,.0f} m3 = {V/crew:,.0f} m3/person; floor {floor:,.0f} m2 = {floor/crew:.0f} m2/person + farms")
    for kk,v in H.items(): print(f"  {kk:50s}{v:12,.0f} t {100*v*1e3/M_hab:5.1f}%")
    print(f"  {'HABITAT TOTAL':50s}{M_hab/1e3:12,.0f} t")
    print("\n=== FUEL PHYSICS (catalysed D-D) ===")
    print(f"energy {E_per_kg:.3e} J/kg ({eps*100:.3f}% of mass); neutrons carry {f_neutron*100:.0f}%")
    print(f"max possible exhaust (100% to jet) {ve_max/1e3:,.0f} km/s; 8,000 km/s needs {eta_needed*100:.1f}% of fusion energy in the jet")
    for nm,e in ENG.items():
        q=2*0.9*sigma*e['T_rad']**4
        print(f"engine case {nm}: absorbs {e['f_abs']*100:.0f}% of fusion power, radiators {e['rad_kgm2']} kg/m2 at {e['T_rad']} K -> effective {sp_eff(**e):.2e} W/kg")
    cases=[(T,'central') for T in (90,100,110,125,150)]+[(100,'optimistic'),(100,'pessimistic'),(125,'pessimistic')]
    for T,nm in cases:
        e=ENG[nm]; sp=sp_eff(**e); r=optimise(T,sp=sp)
        if r is None: print(f"\n=== trip <= {T} yr, engine {nm}: NO SOLUTION"); continue
        Q=e['f_abs']*r['P_fus']; q=2*0.9*sigma*e['T_rad']**4; Arad=Q/q
        print(f"\n=== BEST MISSION, trip <= {T} yr, engine {nm} ({sp:.1e} W/kg) ===")
        print(f" heat absorbed {Q:.2e} W -> radiators {Arad/1e6:,.0f} km2 ({np.sqrt(Arad)/1e3:.0f} km square), {Arad*e['rad_kgm2']/1e9:.3f} Mt")
        print(f" cruise {r['vc']/c:.3f} c ({r['vc']/1e3:,.0f} km/s); start accel {r['a0']:.3f} m/s2")
        print(f" accelerate {r['t1']:.1f} yr over {r['s1']:.2f} ly | cruise {r['t_cr']:.1f} yr | brake {r['t2']:.1f} yr over {r['s2']:.2f} ly | TOTAL {r['T']:.1f} yr")
        print(f" mass: departure {r['M_dep']/1e9:.1f} Mt | fuel {r['fuel']/1e9:.1f} Mt | tanks {r['tanks']/1e9:.2f} Mt | engine {r['M_eng']/1e9:.3f} Mt | at cruise {r['M_cruise']/1e9:.2f} Mt | arrival (dry) {r['M_dry']/1e9:.3f} Mt")
        print(f" thrust {r['F']:.2e} N; jet power {r['P_jet']:.2e} W; fusion power {r['P_fus']:.2e} W; end-of-braking decel {r['a_brake_end']:.3f} m/s2")
        print(f" water to process for deuterium: {r['water']/1e12:,.0f} billion t (an ice body ~{(r['water']/917)**(1/3)/1e3:.1f} km across)")

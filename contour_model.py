"""Unified magsail drag: the shielded region is where the loop's own (vacuum) field in its plane
exceeds the pressure-balance field B_mp = sqrt(mu0 rho v^2). This is the same convention as
Fujita's standoff (far-field dipole gives exactly L = (mu0 I^2 R^4/(16 rho v^2))^(1/6)) and reduces to
the 'doughnut' tube (width 2d, d = mu0 I/(2 pi B_mp)) when d << R. It is continuous in between.
Drag F = Cd (rho v^2/2) * pi (r_out^2 - r_in^2), Cd from Fujita (2004) with x = d_i / Lc,
Lc = r_out (no hole) or (r_out - r_in)/2 (plasma passes through the centre).
"""
import numpy as np
from scipy.special import ellipk, ellipe
mu0=4e-7*np.pi; mp=1.6726e-27

def b_plane(s):
    """B_z in the loop plane at radius s*R, in units of mu0*I/R (signed)."""
    s=np.asarray(s,float); m=4*s/(1+s)**2
    return (1/(2*np.pi))*(1/(1+s))*(ellipk(m)+(1-s**2)/((1-s)**2)*ellipe(m))

# tabulate inverse functions
s_in=np.linspace(0,0.999999,200001); b_in=b_plane(s_in)          # increasing from 0.5 to +inf
s_out=np.geomspace(1.000001,1e4,200001); b_out=np.abs(b_plane(s_out))  # decreasing from inf to ~0
def r_in_of(beta):   # beta = B_mp R/(mu0 I); hole exists if beta > 0.5
    return np.where(beta>0.5, np.interp(beta,b_in,s_in), 0.0)
def r_out_of(beta):
    return np.interp(-np.log(beta), -np.log(b_out), s_out)      # b_out decreasing -> use -log

import os
CDMODEL=os.environ.get('CDMODEL','fujita')   # 'fujita' = Fujita (2004) Eq. 12 (central); 'ashida' = 2.9 refit used in earlier drafts
CDSCALE=float(os.environ.get('CDSCALE','1.0'))  # multiplies Cd (e.g. 4.3/3.6 for Fujita's upper bound)
def cd_fujita(x):
    """Drag coefficient on pi*r_c^2 vs R_L = (ion Larmor radius at r_c)/r_c.
    Fujita (2004) Eq.12: 3.6 exp(-0.28 R^2) for R<1; (3.4/R) exp(-0.22/R^2) for R>=1."""
    x=np.asarray(x,float); xs=np.maximum(x,1e-30)
    if CDMODEL=='ashida':
        cd=np.where(x<1, 2.9*np.exp(-0.35*x**2), 2.4/xs*np.exp(-0.17/xs**2))
    else:
        cd=np.where(x<1, 3.6*np.exp(-0.28*x**2), 3.4/xs*np.exp(-0.22/xs**2))
    return CDSCALE*cd

def drag(v,R,I,n_cm3,coef='fujita'):
    rho=n_cm3*1e6*mp*1.6; di=2.28e5/np.sqrt(n_cm3)
    v=np.asarray(v,float)
    Bmp=np.sqrt(mu0*rho*v**2); beta=Bmp*R/(mu0*I)
    so=r_out_of(beta); si=r_in_of(beta)
    area=np.pi*R**2*(so**2-si**2)
    hole=si>0
    Lc=np.where(hole, R*(so-si)/2, R*so)
    cd=cd_fujita(di/Lc)
    if coef=='freeland': cd=cd*0.345/0.575
    return cd*0.5*rho*v**2*area, ~hole, di/Lc

if __name__=='__main__':
    # checks: far-field dipole limit and thin-tube limit
    I=1e7; R=1e5; n=0.06; v=0.05*2.998e8; rho=n*1e6*mp*1.6
    L=(mu0*I**2*R**4/(16*rho*v**2))**(1/6)
    F,dip,x=drag(v,R,I,n); print("strong-field: contour r_out/L =", np.sqrt(F/(cd_fujita(x)*0.5*rho*v**2*np.pi))/L)
    I=1e6; R=3e7; d=mu0*I/(2*np.pi*np.sqrt(mu0*rho*v**2))
    F,dip,x=drag(v,R,I,n); A=F/(cd_fujita(x)*0.5*rho*v**2); print("thin-tube: area/(4 pi R d) =", A/(4*np.pi*R*d), " d/R =",d/R)
    # continuity scan across the transition
    for I in (1e6,3e6,1e7):
        Rt=(mu0*I**2/(16*rho*v**2))**0.5
        Rs=Rt*np.array([0.5,0.9,0.99,1.01,1.1,2]); F=np.array([drag(v,r,I,n)[0] for r in Rs])
        print(f"I={I/1e6:.0f} MA force vs R/R_switch {np.round(Rs/Rt,2)}: {np.array2string(F,precision=2)}")

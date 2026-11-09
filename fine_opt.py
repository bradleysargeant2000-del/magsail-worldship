"""Sail-assisted optimisation using the SAME fine vc/a0 grid as fusion-only, with cumulative drag integrals.
Includes vf = vc (sail does nothing) so a vanishing sail reproduces fusion-only exactly."""
import numpy as np, warnings; warnings.filterwarnings('ignore')
import full_ship as fs
from contour_model import drag
from scale_model import sail_mass
c=fs.c; yr=fs.yr; D=fs.D_target; Mh=fs.M_hab; dock=fs.dv_dock
sp=fs.sp_eff(**fs.ENG['central']); ve=fs.ve; k=fs.k_tank; ve_eff=ve/(1+k)
caps=(100,125,150)
vcs=np.linspace(0.02,0.07,101)*c          # step 0.0005c (fusion-only optimiser uses 0.0005c too)
a0s=np.geomspace(0.003,0.3,80)
def search(J=3e9,sig=2.9e9,rho_s=1560,Rs=np.geomspace(1e5,1e8,31),Is=np.geomspace(1e5,2e7,31),n=0.06,include_null=True):
    best={T:(np.inf,None) for T in caps}
    vgrid=np.geomspace(1e5,0.075*c,3000)
    for R in Rs:
      for I in Is:
        F,_,_=drag(vgrid,R,I,n)
        G1=np.concatenate([[0],np.cumsum(np.diff(vgrid)*0.5*(1/F[1:]+1/F[:-1]))])
        G2=np.concatenate([[0],np.cumsum(np.diff(vgrid)*0.5*(vgrid[1:]/F[1:]+vgrid[:-1]/F[:-1]))])
        for vc in vcs:
          Fmax=np.interp(vc,vgrid,F); mw,ms,mt=sail_mass(R,I,J=J,sig=sig,rho_s=rho_s,Fmax=Fmax); Ms=mw+ms+mt
          vfs=np.concatenate([np.geomspace(2e5,0.98*vc,25),[vc]])
          for vf in vfs:
            g1=np.interp(vc,vgrid,G1)-np.interp(vf,vgrid,G1); g2=np.interp(vc,vgrid,G2)-np.interp(vf,vgrid,G2)
            Rm=np.exp((vc+vf+dock)/ve_eff)
            den=1-a0s*ve*Rm/(2*sp); ok=den>0
            a=a0s[ok]; Mdep=(Mh+Ms)*Rm/den[ok]; Fe=a*Mdep; Mdry=Mh+Ms+0.5*Fe*ve/sp
            Mbr=Mdry*np.exp((vf+dock)/ve_eff)
            tb=Mbr*g1; sb=Mbr*g2
            mdot=Fe/ve*(1+k); t1=(Mdep-Mbr)/mdot; t2=(Mbr-Mdry)/mdot
            # accel distance (exact): s = ve_eff*( t - (M1/mdot)*ln(M0/M1) ) ... integrate numerically via formula
            s1=ve_eff*(t1 - (Mbr/mdot)*np.log(Mdep/Mbr))
            # braking burn from vf to ~0: s2 = vf*t2 - ve_eff*(t2 - (Mdry/mdot)*ln(Mbr/Mdry))
            s2=vf*t2 - ve_eff*(t2-(Mdry/mdot)*np.log(Mbr/Mdry))
            cr=D-s1-sb-s2; good=cr>0
            T=(t1+np.where(good,cr,0)/vc+tb+t2)/yr
            for cap in caps:
              m=good&(T<=cap)
              if m.any():
                j=np.argmin(np.where(m,Mdep,np.inf))
                if Mdep[j]<best[cap][0]: best[cap]=(Mdep[j],dict(R_km=round(R/1e3),I_MA=round(I/1e6,2),sail_Mt=round(Ms/1e9,4),vc=round(vc/c,4),vf_kms=round(vf/1e3),a0=round(a[j],4),sail_brake_yr=round(tb[j]/yr,1)))
    return best
if __name__=='__main__':
    import sys
    fo={T:fs.optimise(T,sp=sp)['M_dep']/1e9 for T in caps}
    print('fusion only:',{T:round(v,2) for T,v in fo.items()})
    b=search()
    for T in caps:
        print(T,'best with sail:',round(b[T][0]/1e9,2),'Mt',b[T][1])

import sys, numpy as np, fine_opt as fo
f=float(sys.argv[1]); J=float(sys.argv[2]); sig=float(sys.argv[3]); rs=float(sys.argv[4]); tag=sys.argv[5]
# monkeypatch vf choices: sail must remove at least fraction f of cruise speed
orig=np.concatenate
def search(J,sig,rho_s):
    import fine_opt as m
    best={T:(np.inf,None) for T in m.caps}
    vgrid=np.geomspace(1e5,0.075*m.c,3000)
    import os
    for R in np.geomspace(float(os.environ.get('RLO','1e5')),float(os.environ.get('RHI','1e8')),int(os.environ.get('NR','31'))):
      for I in np.geomspace(float(os.environ.get('ILO','1e5')),float(os.environ.get('IHI','2e7')),int(os.environ.get('NI','31'))):
        import os; F,_,_=m.drag(vgrid,R,I,float(os.environ.get("NDEN","0.06")))
        G1=np.concatenate([[0],np.cumsum(np.diff(vgrid)*0.5*(1/F[1:]+1/F[:-1]))])
        G2=np.concatenate([[0],np.cumsum(np.diff(vgrid)*0.5*(vgrid[1:]/F[1:]+vgrid[:-1]/F[:-1]))])
        for vc in m.vcs:
          Fmax=np.interp(vc,vgrid,F); mw,ms,mt=m.sail_mass(R,I,J=J,sig=sig,rho_s=rho_s,Fmax=Fmax); Ms=mw+ms+mt
          for vf in np.geomspace(2e5,(1-f)*vc,20):
            g1=np.interp(vc,vgrid,G1)-np.interp(vf,vgrid,G1); g2=np.interp(vc,vgrid,G2)-np.interp(vf,vgrid,G2)
            Rm=np.exp((vc+vf+m.dock)/m.ve_eff); den=1-m.a0s*m.ve*Rm/(2*m.sp); ok=den>0
            a=m.a0s[ok]; Mdep=(m.Mh+Ms)*Rm/den[ok]; Fe=a*Mdep; Mdry=m.Mh+Ms+0.5*Fe*m.ve/m.sp
            Mbr=Mdry*np.exp((vf+m.dock)/m.ve_eff); tb=Mbr*g1; sb=Mbr*g2
            mdot=Fe/m.ve*(1+m.k); t1=(Mdep-Mbr)/mdot; t2=(Mbr-Mdry)/mdot
            s1=m.ve_eff*(t1-(Mbr/mdot)*np.log(Mdep/Mbr)); s2=vf*t2-m.ve_eff*(t2-(Mdry/mdot)*np.log(Mbr/Mdry))
            cr=m.D-s1-sb-s2; good=cr>0; T=(t1+np.where(good,cr,0)/vc+tb+t2)/m.yr
            for cap in m.caps:
              mm=good&(T<=cap)
              if mm.any():
                j=np.argmin(np.where(mm,Mdep,np.inf))
                if Mdep[j]<best[cap][0]: best[cap]=(Mdep[j],(round(R/1e3),round(I/1e6,2),round(Ms/1e9,3),round(vc/m.c,4),round(vf/1e3),round(tb[j]/m.yr,1)))
    return best
b=search(J,sig,rs)
print(tag,'f=',f,{T:(round(v[0]/1e9,1) if v[1] else None, v[1]) for T,v in b.items()},flush=True)

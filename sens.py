from full_ship import *
import full_ship as fs
def line(tag,r):
    if r is None: print(f"{tag:55s} no solution"); return
    print(f"{tag:55s} depart {r['M_dep']/1e9:7.1f} Mt | cruise {r['vc']/c:.3f}c | dry {r['M_dry']/1e9:.2f} Mt | accel {r['t1']:.1f} yr | fuel {r['fuel']/1e9:.1f} Mt")
spc=sp_eff(**ENG['central']); spp=sp_eff(**ENG['pessimistic'])
for T in (150,200,300): line(f"pessimistic engine, trip<={T}",optimise(T,sp=spp))
for v in (5e6,1e7,1.2e7):
  for T in (100,125): line(f"central engine, ve {v/1e3:.0f} km/s, trip<={T}",optimise(T,sp=spc,ve=v))
# heavier shield (7 t/m2 meets 20 mSv/yr; central case uses 6.5)
Mh2=M_hab+H['Radiation shield (water 6.5 t/m2)']*1e3*(7/6.5-1)
for T in (100,125): line(f"central, shield 7 t/m2 (hab {Mh2/1e9:.2f} Mt), trip<={T}",optimise(T,sp=spc,M_hab=Mh2))

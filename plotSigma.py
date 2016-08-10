from fermion2PointBN import A, zReal, z, makeCache, dispersion, G, crossing
import pylab as pl
import numpy as np
from numpy import pi, sqrt, arccosh, sign

N=1500
v=.5

kxs=[0,.02,-.1,.1,crossing(v)]
omegas=np.linspace(-1, 1, N)

def G1loop(omega, kx, v):
    return -((4 *pi*sqrt(1 - v**2))/(4*pi*sqrt(1 - v**2)* (kx* v - omega) - (arccosh(0j+( kx - v *omega)/(kx*v - omega)) if kx>0 else -(arccosh(0j+( kx - v *omega)/(kx*v - omega))).conjugate())))

for ki in range(len(kxs)):
    kx=kxs[ki]
    kxt=['0','0.02','-.1','0.1','k_x^*'][ki]
    Gs=[G(omega, kx, v) for omega in omegas]
    G1s=[G1loop(omega, kx, v) for omega in omegas]
    import matplotlib.pyplot as plt
    from fig import fig, saveFig, fill_between
    for i in [0,1]:
        fig(2*ki+i,size=8.25)
        pl.grid()
        pl.plot(omegas,[[g.real,-2*g.imag][i] for g in Gs],label='$'+['\mathrm{Re}(G_R)','A'][i]+'$',ls='-',zorder=3,color='k')
        pl.plot(omegas,[[g.real,-2*g.imag][i] for g in G1s],label=r'$'+['\mathrm{Re}(G_R^{\mathrm{1loop}})','A^{\mathrm{1loop}}'][i]+'$',ls='--',zorder=3,color='r')
        if i and ki:
            omegaPole=omegas[G1s.index(max(G1s))]-.5*(omegas[1]-omegas[0])
            pl.plot([omegaPole, omegaPole], [0,1000],ls='--',lw=1.5,zorder=3,color='r')
        pl.ylim([-5,40] if i else [-20,20])
        pl.legend(loc=2)
        pl.xlabel(r'$\omega/\lambda^2$')
        pl.ylabel(r'$\lambda^2'+['G','A'][i]+'(\omega, k_x)$')
        saveFig(['re','im'][i]+'G_v='+str(v)+',kx='+str(kx))
plt.show()

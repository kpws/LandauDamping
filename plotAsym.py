from results import loadG, getGfun, Gquenchedv1, GlargeNf
import sys
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from fig import fig, saveFig, fill_between, grid_selfmade

Nfs=[.01,1,100]
L=2000
name='run10'
ws=[.01,1,20]

Gobs=[loadG(name, Nf, L, 15) for Nf in Nfs]
Gs=[getGfun(Gob, 'cubic') for Gob in Gobs]
wm=Gobs[0][0]
n=len(Gobs[0][1])
ws=np.linspace(-wm,wm,n)[n/2:]
skip=1
ws=np.concatenate((ws[skip:n/256],np.logspace(np.log10(ws[n/256]),np.log10(ws[-1]),20)))
print('Datapoints: '+str(n))
N=2000
wmin=wm/(n-1)
wmax=wm
xs=wmin*np.exp(np.linspace(0, np.log(wmax/wmin), N))

colors=['red', 'green', 'blue', 'violet']

f=fig(2,size=14)
phis=np.linspace(0,2*np.pi,20)[:-1]
rs=[.05,1]
j=0
for i in range(len(rs)):
	phases=np.exp(1j*phis)*rs[i]
	pl.plot(phis, [np.angle(Gs[j](p.imag,p.real)) for p in phases])
pl.show()
sys.exit(0)


f=fig(0,size=14)
pl.plot(xs, np.abs([Gquenchedv1(0,w) for w in xs]),label=(r'$N_fk_f=0$'),color='black')
pl.plot(xs, np.abs([1/w for w in xs]),label=(r'$\lambda=0$'),color='black',ls='--')
for j in range(len(Nfs)):
	pl.loglog(ws, abs(Gs[j](0,ws)),'+',label='$N_fk_f='+str(Nfs[j])+r'\lambda^2$',color=colors[j])
#for j in range(len(Nfs)):
#	pl.plot(xs, np.abs([GlargeNf(Nfs[j], 0, w) for w in xs]),label=(r'$N_f\rightarrow\infty$'),color=colors[j],ls='-')
#, N_fk_f='+str(Nfs[j])+'\lambda^2
		#if r==1:
		#	pl.plot(xs, cf(SlargeNf(Nfs[j],xs)),color=colors[j],linestyle='--')
#pl.plot([],[],label=r'$k_x/\lambda^2='+wss+'$',color='white')
pl.xlim([wmin,wmax])
pl.ylim([7e-2,2e2])

grid_selfmade(f.get_axes()[0])

pl.legend(loc=1)
pl.xlabel(r'$k_x/\lambda^2$')
pl.ylabel(r'$|G(\omega=0, k_x)|/\lambda^2$')

saveFig('absG_vs_ks_L='+str(L))

#pl.show()
#sys.exit(0)

f=fig(1,size=14)
pl.plot(xs, np.abs([Gquenchedv1(w,0) for w in xs]),label=(r'$N_fk_f=0$'),color='black')
pl.plot(xs, np.abs([1/w for w in xs]),label=(r'$\lambda=0$'),color='black',ls='--')
for j in range(len(Nfs)):
	pl.loglog(ws, abs(Gs[j](ws,0)),'+',label='$N_fk_f='+str(Nfs[j])+r'\lambda^2$',color=colors[j])
for j in range(len(Nfs)):
	pl.plot(xs, np.abs([GlargeNf(Nfs[j], w, 0) for w in xs]),label=(r'$N_f\rightarrow\infty$'),color=colors[j],ls='-')
#, N_fk_f='+str(Nfs[j])+'\lambda^2
		#if r==1:
		#	pl.plot(xs, cf(SlargeNf(Nfs[j],xs)),color=colors[j],linestyle='--')
#pl.plot([],[],label=r'$k_x/\lambda^2='+wss+'$',color='white')
pl.xlim([wmin,wmax])
pl.ylim([7e-2,2e2])

grid_selfmade(f.get_axes()[0])

pl.legend(loc=1)
pl.xlabel(r'$\omega/\lambda^2$')
pl.ylabel(r'$|G(\omega, k_x=0)|/\lambda^2$')

saveFig('absG_vs_omega_L='+str(L))

pl.show()

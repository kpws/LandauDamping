from results import loadG, loadGBenchmark, getSfun, getGfun
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from fig import fig, saveFig, fill_between, grid_selfmade

Nfs=[.01,.1,1]
#ws=[.05,1,12]

#L=8000
#name='4'
#ws=[.01,.1,3]

L=2000
name='benchmark1'
ws=[.1,1,10]

Gobs=[loadG(name, Nf, L, 15) for Nf in Nfs]
Ss=[getSfun(Gob) for Gob in Gobs]
print('Datapoints: '+str(len(Gobs[0][1])))

nameR='benchmarkRef3'
GobRs=[loadGBenchmark(nameR, Nf) for Nf in Nfs]
SsR=[getSfun(GobR) for GobR in GobRs]
print('DatapointsR: '+str(len(GobRs[0][1])))

N=1500
wm=min(Gobs[0][0], GobRs[0][0])
xs=np.linspace(-wm, wm, N)
wss=reduce(lambda a,b:a+', '+b,map(str,ws))

colors=['red', 'green', 'blue', 'violet']


print Ss[0](1,.1)
print SsR[0](1,.1)

for r in [0,1]:
	f=fig(r,size=14)
	for i in range(len(ws)):
		kx=ws[i]
		cf=[np.real,np.imag][r]
		for j in range(len(Nfs)):
			pl.plot(xs, cf(Ss[j](xs,kx)),label='$N_fk_f='+str(Nfs[j])+r'\lambda^2$' if i==0 else '',color='black')
			pl.plot(xs, cf(SsR[j](xs,kx)),label='$N_fk_f='+str(Nfs[j])+r'\lambda^2$' if i==0 else '',color=colors[j])
	pl.xlim([-wm,wm])
	if r==1:
		pl.ylim([-.11, .11])
	else:
		pl.ylim([-0,.09])
	grid_selfmade(f.get_axes()[0])

	pl.legend(loc=4)
	pl.xlabel(r'$\omega/\lambda^2$')
	pl.ylabel(r'$\mathrm{'+['Re','Im'][r]+'}(\Sigma(\omega, k_x))/\lambda^2$')

	saveFig(['re','im'][r]+'_Sigma_BM_vs_omega')

plt.show()
for r in [0,1]:
	f=fig(2+r,size=14)
	for i in range(len(ws)):
		w=ws[i]
		cf=[np.real,np.imag][r]
		for j in range(len(Nfs)):
			pl.plot(xs, cf(Ss[j](w,xs)),label='$N_fk_f='+str(Nfs[j])+r'\lambda^2$' if i==0 else '',color=colors[j])
			#iif r==1:
			#	pl.plot([-wm,wm], 2*[cf(SlargeNf(Nfs[j],w))],color=colors[j],linestyle='--')
		pl.plot(xs, cf([Squenchedv1(w,kx) for kx in xs]),label=(r'$N_fk_f=0$' if i==0 else ''),color='black')
	pl.plot([],[],label=r'$\omega/\lambda^2='+wss+'$',color='white')
	pl.xlim([-wm,wm])
	if r==0:
		pl.ylim([-.11, .11])
	else:
		pl.ylim([-0,.09])
	grid_selfmade(f.get_axes()[0])

	pl.legend(loc=4)
	pl.xlabel(r'$k_x/\lambda^2$')
	pl.ylabel(r'$\mathrm{'+['Re','Im'][r]+'}(\Sigma(\omega, k_x))/\lambda^2$')

	saveFig(['re','im'][r]+'_Sigma_vs_kx_L='+str(L))

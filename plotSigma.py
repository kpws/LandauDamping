from results import loadG, getSfun, Squenchedv1, SlargeNf
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from fig import fig, saveFig, fill_between, grid_selfmade
'''
Nfs=[.01,.1,1,10]
L=2000
name='run10'
ws=[.01,1,20]
n2=15
'''
Nfs=[.01]
L=10000
name='run_sub1'
ws=[.01,1,20]
aa=True
n2=11

Gobs=[loadG(name, Nf, L, n2) for Nf in Nfs]
Ss=[getSfun(Gob, addAsym=aa) for Gob in Gobs]
wm=Gobs[0][0]
print('Datapoints: '+str(len(Gobs[0][1])))
N=200
xs=np.linspace(-wm, wm, N)
wss=reduce(lambda a,b:a+', '+b,map(str,ws))

colors=['red', 'green', 'blue', 'violet']

for r in [0,1]:
	f=fig(r,size=14)
	for i in range(len(ws)):
		kx=ws[i]
		cf=[np.real,np.imag][r]
		for j in range(len(Nfs)):
			pl.plot(xs, [cf(Ss[j](x,kx)) for x in xs],label='$N_fk_f='+str(Nfs[j])+r'\lambda^2$' if i==0 else '',color=colors[j])
			if r==1:
				pl.plot(xs, cf(SlargeNf(Nfs[j],xs)),color=colors[j],linestyle='--')
		pl.plot(xs, cf([Squenchedv1(w,kx) for w in xs]),label=(r'$N_fk_f=0$' if i==0 else ''),color='black')
	pl.plot([],[],label=r'$k_x/\lambda^2='+wss+'$',color='white')
	pl.xlim([-wm,wm])
	if r==1:
		pl.ylim([-.11, .11])
	else:
		pl.ylim([-0,.09])
	grid_selfmade(f.get_axes()[0])

	pl.legend(loc=4)
	pl.xlabel(r'$\omega/\lambda^2$')
	pl.ylabel(r'$\mathrm{'+['Re','Im'][r]+'}(\Sigma(\omega, k_x))/\lambda^2$')

	saveFig(['re','im'][r]+'_Sigma_vs_omega_L='+str(L))

for r in [0,1]:
	f=fig(2+r,size=14)
	for i in range(len(ws)):
		w=ws[i]
		cf=[np.real,np.imag][r]
		for j in range(len(Nfs)):
			pl.plot(xs, [cf(Ss[j](w,x)) for x in xs],label='$N_fk_f='+str(Nfs[j])+r'\lambda^2$' if i==0 else '',color=colors[j])
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
plt.show()

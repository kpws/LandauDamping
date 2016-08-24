from results import loadG, loadGBenchmark, getSfun, getGfun, GlargeNfApp
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from fig import fig, saveFig, fill_between, grid_selfmade

Nfs=[.01,.1,1,10,100]
L=250
name='benchmark7'
ws=[1]
n2=13

Nfs=[.01,.1,1,10,100]
L=500
name='benchmark9'
ws=[1]
n2=14

Nfs=[.1]
L=500
name='benchmark_sub1'
ws=[1]
n2=12


Gobs=[loadG(name, Nf, L, n2) for Nf in Nfs]
Ss=[getSfun(Gob, addFree=True) for Gob in Gobs]
nameR='benchmarkRef3'
#GobRs=[loadGBenchmark(nameR, Nf) for Nf in Nfs]
SsR=[lambda w, kxm, Nf=Nf: -1/GlargeNfApp(Nf, w, kx)+1j*w-kx for Nf in Nfs]
#SsR=[getSfun(GobR) for GobR in GobRs]
#print('DatapointsR: '+str(len(GobRs[0][1])))
#print('Max omegaR: '+str(GobRs[0][0]))

N=2000
wm=Gobs[0][0]
print('Datapoints: '+str(len(Gobs[0][1])))
print('Max omega: '+str(Gobs[0][0]))
xs=np.linspace(-wm, wm, N)
Nsparse=60
xss=np.linspace(-wm, wm, Nsparse)
wss=reduce(lambda a,b:a+', '+b,map(str,ws))

colors=['red', 'green', 'blue', 'violet','yellow']


#print Ss[0](1,.3)
#print SsR[0](1,.3)

for r in [0,1]:
	f=fig(r,size=12)
	for i in range(len(ws)):
		kx=ws[i]
		cf=[lambda x:-np.real(x),np.imag][r]
		for j in range(len(Nfs)):
			pl.plot(xs, cf(SsR[j](xs,kx)),color='black')
			pl.plot(xss, cf(Ss[j](xss,kx)),'+',label='$N_fk_f='+('' if Nfs[j]==1 else str(Nfs[j]))+r'\lambda^2$' if i==0 else '',color=colors[j])
	pl.xlim([-wm,wm])
	pl.plot([],[],label=r'$\mathrm{Analytic\ solutions}$',color='black')
	if r==1:
		pl.ylim([-.21, .21])
		pl.legend(loc=1)
	else:
		pl.ylim([0,.0028])
		pl.legend(loc=1)
	grid_selfmade(f.get_axes()[0])

	pl.xlabel(r'$\omega/\lambda^2$')
	pl.ylabel(r'$\mathrm{'+['-Re','Im'][r]+'}(\Sigma(\omega, k_x))/\lambda^2$')

	saveFig(['re','im'][r]+'_Sigma_BM_vs_omega')

plt.show()

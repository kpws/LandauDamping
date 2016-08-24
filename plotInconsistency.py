from results import loadG, loadGBenchmark, getSfun, getGfun, GlargeNfApp
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from fig import fig, saveFig, fill_between, grid_selfmade

N=35000
wm=.021234;

colors=['red', 'green', 'blue', 'violet','yellow']


f=fig(0,size=10)
Nf=1.
kx=1
xs=np.linspace(kx-wm, kx+wm, N)
pl.plot(xs, -2*np.imag(GlargeNfApp(Nf,-1j*xs+0.000001,kx)),label='',color='black')
pl.plot(np.linspace(kx-.0005,kx+.0005,1000),500*[-120000,120000],color='black')
#pl.plot(xs, -np.imag(GlargeNfApp(Nf,-1j*xs+0.00001,kx)),label='',color='blue')
#pl.yscale('log')
pl.xlim([xs[0],xs[-1]])
#pl.legend(loc=4)
pl.ylim([-100000,100000])
grid_selfmade(f.get_axes()[0])
pl.xlabel(r'$\omega/\lambda^2$')
pl.ylabel(r'$A(\omega, k_x=\lambda^2)/\lambda^2$')

saveFig('inconsistency')

plt.show()

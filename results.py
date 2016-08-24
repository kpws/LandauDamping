from numpy import sqrt, pi, linspace
import numpy as np
from scipy import interpolate

def Squenchedv1(w,kx):
	if 0<w.real:
		sgn=1
	else:
		sgn=-1
		w=-w.conjugate()

	kxA=abs(kx)
	A=-8*kxA**3*pi**2 + 3j*w + 24j*kxA**2*pi**2*w - 8j*pi**2*w**3 + 3*kxA*(1 + 8*pi**2*w**2) + sqrt(9*kxA**2 - 48*kxA**4*pi**2 + 96j*kxA**3*pi**2*w - 9*w**2 + 48*pi**2*w**4 + 6j*kxA*w*(3 + 16*pi**2*w**2))
	S=(4*pi**(10./3)*(kxA - 1j*w) + (4*pi**4*(kxA - 1j*w)**2*(1 - 1j*sqrt(3)))/A**(1./3) + A**(1./3)*pi**(8./3)*(1 + 1j*sqrt(3)))/(8.*pi**(10./3))
	if kx<0: S=-S.conjugate()
	if sgn<0: S=S.conjugate()
	return S

from scipy import special
from mpmath import hyp1f2
Gamma2o3=special.gamma(2./3)
Gamma1o3=special.gamma(1./3)
Gammam1o3=special.gamma(-1./3)

def GlargeNfApp(Nf, w, kx, subLeading=False):
	sgn = np.sign(w.real)
	w=w*sgn
	kx=kx*sgn
	l0=(3.*np.sqrt(3.))**3*(2*np.pi)**5*Nf/(2*Gamma2o3)**3
	p1 = np.cos(w/(np.sqrt(l0)*(w+1j*kx)**(3./2)))/(1j*w-kx)
	arg = - w**2/(4*l0*(w+1j*kx)**3)
	#print('w: '+str(w))
	#print('arg: '+str(arg))
	p2 = 6j*np.sqrt(3)*Gamma1o3*w**(2./3)/(8*np.pi*l0**(1./3)*(w+1j*kx)**2) * complex(hyp1f2(1.,5./6,4./3, arg))  
	p3 = 3j*np.sqrt(3)*Gammam1o3*w**(4./3)/(8*np.pi*l0**(2./3)*(w+1j*kx)**3) * complex(hyp1f2(1., 7./6, 5./3, arg)) 
	if subLeading: sgn=sgn*np.exp(3/(36*np.pi*np.pi)/np.sqrt(Nf))
	return sgn*(p1+p2+p3)

def oldGlargeNfApp(Nf, w, kx, subLeading=False):
	sgn = np.sign(w.real)
	w=w*sgn
	kx=kx*sgn
	l0=(3.*np.sqrt(3.))**3*(2*np.pi)**5*Nf/(2*Gamma2o3)**3
	p1 = np.cos(w/(np.sqrt(l0)*(w+1j*kx)**(3./2)))/(1j*w-kx)
	arg = - w**2/(4*l0*(w+1j*kx)**3)
	print arg
	p2 = 6j*np.sqrt(3)*Gamma1o3*w**(2./3)/(8*np.pi*l0**(1./3)*(w+1j*kx)**2) * np.array([complex(hyp1f2(1.,5./6,4./3, a))  for a in arg]) 
	p3 = 3j*np.sqrt(3)*Gammam1o3*w**(4./3)/(8*np.pi*l0**(2./3)*(w+1j*kx)**3) * np.array([complex(hyp1f2(1., 7./6, 5./3, a)) for a in arg])
	if subLeading: sgn=sgn*np.exp(3/(36*M_PI*M_PI)/np.sqrt(Nf))
	return sgn*(p1+p2+p3)

def Gquenchedv1(w,kx):
	return 1/(1j*w - kx + Squenchedv1(w,kx))

def SlargeNf(Nf, w):
	return 1j*np.sign(w)*np.absolute(w)**(2./3)/((2*pi)**(5./3)*sqrt(3)*Nf**(1./3))

def G0(w, kx):
	return 1/(1j*w - kx)

def GlargeNf(Nf, w, kx):
	return 1/(1j*w - kx + SlargeNf(Nf,w))

def loadGBenchmark(name, Nf):
	with open('benchmarkRefOutput/'+name+'_Nf='+str(Nf), 'rb') as f:
		n = np.fromfile(f, dtype=np.int32, count=1)[0]
		wm = np.fromfile(f, dtype=np.float64, count=1)[0]
		v = np.fromfile(f, dtype=np.complex128, count=n*n).reshape((n,n)).transpose()
	return (wm,v)

def loadG(name, Nf, L, n2):
	n=2**n2;
	with open('cgsOutput/cgsOutput_'+name+'_'+str(Nf)+'_'+str(n)+'_'+str(L), 'rb') as f:
		t=np.int64
		Nf = np.fromfile(f, dtype=np.float64, count=1)[0]
		L = np.fromfile(f, dtype=np.float64, count=1)[0]
		np.fromfile(f, dtype=t, count=1)
		np.fromfile(f, dtype=np.float64, count=1)
		np.fromfile(f, dtype=np.float64, count=1)
		np.fromfile(f, dtype=t, count=1)
		n = np.fromfile(f, dtype=t, count=1)[0]
		wm = np.fromfile(f, dtype=np.float64, count=1)[0]
		v = np.fromfile(f, dtype=np.complex128, count=n*n).reshape((n,n))
	return (wm,v, Nf)

def getGfun(Gob, k='cubic'):
	n=len(Gob[1])
	ws=linspace(-Gob[0],Gob[0],n).reshape((1,n))
	#scipy.interpolate.RectBivariateSpline is perhaps faster
	return lambda w,kx:interpolate.interp2d(ws,ws,Gob[1].real,kind='cubic')(w,kx) + 1j*interpolate.interp2d(ws,ws,Gob[1].imag,kind=k)(w,kx)

def getSfun(Gob, addFree=False, addAsym=False):
	n=len(Gob[1])
	ws=linspace(-Gob[0],Gob[0],n).reshape((1,n))
	iG0v=1j*ws-ws.transpose()
	if addAsym:
		return lambda w,kx: -1./(0*GlargeNfApp(Gob[2], w, kx, subLeading=True) + interpolate.interp2d(ws,ws,Gob[1].real,kind='cubic')(w,kx) + 1j*interpolate.interp2d(ws,ws,Gob[1].imag,kind='cubic')(w,kx))+1./G0(w,kx)
	add=0.
	if addFree: add=1./iG0v
	Sv=-1/(Gob[1]+add) + iG0v
	return lambda w,kx:interpolate.interp2d(ws,ws,Sv.real,kind='cubic')(w,kx) + 1j*interpolate.interp2d(ws,ws,Sv.imag,kind='cubic')(w,kx)


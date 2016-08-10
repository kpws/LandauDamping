from numpy import sqrt, pi

def Squenchedv1(w,kx):
	kxA=abs(kx)
	A=-8*kxA**3*pi**2 + 3*I*w + 24*I*kxA**2*pi**2*w - 8*I*pi**2*w**3 + 3*kxA*(1 + 8*pi**2*w**2) + sqrt(9*kxA**2 - 48*kxA**4*pi**2 + 96*I*kxA**3*pi**2*w - 9*w**2 + 48*pi**2*w**4 + 6*I*kxA*w*(3 + 16*pi**2*w**2))
	S=(4*pi**(10./3)*(kxA - I*w) + (4*pi**4*(kxA - I*w)**2*(1 - I*sqrt(3)))/A**(1./3) + A**(1./3)*pi**(8./3)*(1 + I*sqrt(3)))/(8.*pi**(10./3))
	if kx<0:S=-conj(S)

def Gquenchedv1(w,kx):
	return 1/(I*w - kx + Squenchedv1(w,kx))

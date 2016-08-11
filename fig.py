from pylab import figure, savefig, rcParams, text, gca
import pylab as pl
from math import sqrt, atan, pi
import numpy as np

lineWidth=.75
font=10 #12

params = {'backend': 'ps',
          'axes.labelsize': font,
          'font.size': font,
          'legend.fontsize': font,
          'xtick.labelsize': font,
          'ytick.labelsize': font,
          'text.usetex': True,
          'lines.linewidth':lineWidth,
          #'legend.linewidth':0.5,
          'axes.linewidth':lineWidth,
          'grid.linewidth':.5,
          'patch.linewidth':lineWidth
          }

rcParams.update(params)


def grid_selfmade(ax,minor=False):
	y_axis=ax.yaxis.get_view_interval()
	for xx in ax.xaxis.get_ticklocs():
		ax.plot([xx,xx],y_axis,linestyle=':',color='black',linewidth=0.5,zorder=0)
	if minor==True:
		for xx in ax.xaxis.get_ticklocs(minor=True):
			ax.plot([xx,xx],y_axis,linestyle=':',color='black',linewidth=0.5,zorder=0)
	x_axis=ax.xaxis.get_view_interval()
	for yy in ax.yaxis.get_ticklocs():
		ax.plot(x_axis,[yy,yy],linestyle=':',color='black',linewidth=0.5,zorder=0,)
	if minor==True:
		for yy in ax.yaxis.get_ticklocs(minor=True):
			ax.plot(x_axis,[yy,yy],linestyle=':',color='black',linewidth=0.5,zorder=0)

def fig(i,size=14,aspect=None):
    c=0.393701#inches per cm
    if aspect==None:
        aspect=(1+sqrt(5))/2
    return figure(i,figsize=(size*c,size*c/aspect))

def saveFig(name):
    savefig('figs/'+name+'.pdf', bbox_inches='tight')
'''import plotText
def printText(x,y,m,k,txt):
    a=m+k*x[0]-y[0]
    for i in range(1,len(x)):
        if (m+k*x[i]-y[i])*a<0:
            k2=(y[i]-y[i-1])/(x[i]-x[i-1])
            m2=y[i]-k2*x[i]
            x=(m-m2)/(k2-k)
            y=m+k*x
            rot=atan(k2)/(2*pi)*360
            rot=gca().transData.transform_angles(np.array([rot]),np.array([[x,y]]))[0]
            #text(x,y,txt,rotation=rot,ha='center',va='center',size=8,backgroundcolor='white')
            plotText.text(x,y,rot,txt,size=8)
            return
    print('Never crossed line...')
'''
def fill_between(x, y1, y2=0, ax=None, **kwargs):
    """Plot filled region between `y1` and `y2`.

    This function works exactly the same as matplotlib's fill_between, except
    that it also plots a proxy artist (specifically, a rectangle of 0 size)
    so that it can be added it appears on a legend.

    """
    ax = ax if ax is not None else pl.gca()
    ax.fill_between(x, y1, y2, **kwargs)
    p = pl.Rectangle((0, 0), 0, 0, **kwargs)
    ax.add_patch(p)
    return p

def getPlotY(x0,x1,f,yf,maxTurn=0.1,minN=50,maxN=0):
    xs=np.linspace(x0,x1,minN)
    v=[]
    for x in xs:
        o=f(x)
        v.append((x,yf(o),o))
    cosMaxTurn=np.cos(maxTurn)
    change=True
    while change:
        change=False
        i=0
        while i<len(v)-2:
            cosTurn=(  ((v[i+1][0]-v[i][0])*(v[i+2][0]-v[i+1][0])+ (v[i+1][1]-v[i][1])*(v[i+2][1]-v[i+1][1]))/
                         np.sqrt(((v[i+1][0]-v[i  ][0])**2+(v[i+1][1]-v[i  ][1])**2)*
                          ((v[i+2][0]-v[i+1][0])**2+(v[i+2][1]-v[i+1][1])**2))  )
            if cosTurn<cosMaxTurn and (maxN==0 or maxN>len(v)):
                  #print(np.arccos(cosTurn))
                  mx=(v[i+1][0]*1+v[i+2][0])/2
                  o=f(mx)
                  v.insert(i+2,(mx,yf(o),o))
                  mx=(v[i][0]+v[i+1][0]*1)/2
                  o=f(mx)
                  v.insert(i+1,(mx,yf(o),o))
                  i+=3
                  change=True
            i+=1
    return zip(*v)

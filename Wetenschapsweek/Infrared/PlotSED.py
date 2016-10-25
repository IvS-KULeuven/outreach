import numpy as np
import matplotlib.pyplot as plt, matplotlib as mpl, numpy as np, sys, math
from ivs.io import ascii
from matplotlib.widgets import Slider, Button, RadioButtons
mpl.rc('font',size=17)
mpl.rc('mathtext',fontset='custom')
mpl.rc('mathtext',default='sf')

def radiation(wavelengths,temperature):
    h=6.626068
    k=1.380650e-7 
    c=2.99792458 
    length=len(wavelengths)
    blackbody_rad = (1e10)*4.*math.pi*h*(c**2)*(1000./wavelengths)**5*(np.ones(length)/(np.exp(h*c*np.ones(length)/(k*temperature*wavelengths))-np.ones(length)))
    return blackbody_rad


wa = 10.**(np.linspace(-1., 2.7, 1000))
wa2 = 10.**(np.linspace(-1., 2.7, 1000))
star = sys.argv[1]

if(star == 'BD+48_1220_shell'):
  fct = 2.2*10.**(-20.)
  Tstar = 7900.
  Tmin = 1.
  Tmax = 1000.0
  ratmin = 0.01
  ratmax = 5000.
  Rrat = 100.
  Td0 = 10.
elif(star == 'TW_Cam_disc'):
  fct = 2.2*10.**(-19.)
  Tstar = 4800.
  Tmin = 1.
  Tmax = 1600.0
  ratmin = 0.01
  ratmax = 100.
  Td0 = 10.
else:
  fct = 10.**(-20.)
  Tstar = 5778.
  Tmin = 1.
  Tmax = 1600.0
  ratmin = 0.01
  ratmax = 5000.
  Td0 = 10.

Rrat = 0.5*(Tstar/Td0)**2.
data = ascii.read2recarray(star+'.phot')
lmax = 2.89777*10.**3. / Td0

fig = plt.figure(1,figsize=(8,8))
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.15, bottom=0.25)
plt.title(star)
ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')
plt.errorbar(data['cwave'],data['cwave']*data['cmeas'],yerr=data['cwave']*data['e_cmeas'],fmt='ko')
ax.set_xlabel('Wavelength $\lambda$ [$\AA$]')
ax.set_ylabel('$\lambda f_{\lambda}$ [erg/s/$cm^2]')
plt.plot(wa*10.**4.,radiation(wa*10.**4.,Tstar)*wa*10.**4.*fct,'r--',lw=2.)
plt.annotate('$T_{\sf star}$ = %4iK'%Tstar,xy=(0.54,0.93),xycoords='axes fraction',xytext=(0.54,0.93),textcoords='axes fraction')
d =plt.text(1.45*10.**5.,4.0*10.**(-9.),'$R_{\sf dust}/R_{\sf star}$ = %06.1f'%Rrat)
e =plt.text(1.45*10.**5.,2.5*10.**(-9.),'$\lambda_{\sf max}$ = %4.1f$\mu$m'%lmax)

s = (Rrat**2.)*fct*radiation(wa2*10.**4.,Td0)*wa2*10.**4.
l,=plt.plot(wa2*10.**4.,s,'b-',lw=2.)
k,=plt.plot(wa2*10.**4.,radiation(wa2*10.**4.,Tstar)*wa2*10.**4.*fct+s,'r-',lw=2.)

plt.axis([10.**3., 10.**7., 10.**(-11.), 10.**(-8.)])

axcolor = 'white'
axfreq = plt.axes([0.15, 0.1, 0.65, 0.03], axisbg=axcolor)

sfreq = Slider(axfreq, '$T_{\sf dust}$ [K]', Tmin, Tmax, valinit=Td0)


def update(val):
    Tdust = sfreq.val
    Rratio = 0.5*(Tstar/Tdust)**2.
    lmaxn = 2.89777*10.**3. / Tdust
    l.set_ydata((Rratio**2.)*fct*radiation(wa2*10.**4.,Tdust)*wa2*10.**4.)
    k.set_ydata(radiation(wa2*10.**4.,Tstar)*wa2*10.**4.*fct+(Rratio**2.)*fct*radiation(wa2*10.**4.,Tdust)*wa2*10.**4.)
    d.set_text('$R_{\sf dust}/R_{\sf star}$ = %6.1f'%Rratio)
    e.set_text('$\lambda_{\sf max}$ = %4.1f$\mu$m'%lmaxn)
    fig.canvas.draw_idle()
sfreq.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sfreq.reset()
button.on_clicked(reset)


plt.show()
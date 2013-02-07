import numpy as np
from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons


# initial values and constants

G = 6.67*10**(-11)
H0 = 0.0753
mass2 = 1.*10**(-40)

a0 = 1.

adot0 = H0

phi0 = 1.

def init_phidot(phi, m2):
    return np.sqrt((0.525/(np.pi*G))*H0**2 - m2*phi**2)

phidot0 = init_phidot(phi0, mass2)

dt = 0.00001


# stupid euler routines

def D_phidot(a, adot, f, fdot, m2):
    return fdot + (3*(adot/a)*fdot + m2*f)*dt

def D_phi(f, fdot):
    return f - fdot*dt

def D_a(a, adot):
    return a - adot*dt

def Friedman(a, f, fdot, m2):
    return np.sqrt((0.3*H0**2)/a + ((4./3.)*np.pi*G)*(fdot**2 + m2*f**2)*a**2)

# have at it

def run_dat_shiet(a, adot, phi, phidot, m2):
    k = 1
    alist = [a]
    adotl = [adot]
    phil = [phi]
    phidotl = [phidot]
    t = 0.
    tlist = [0.]

    while a > 0.0001 and k < 10000000:

        fdot_old = phidot
        phidot = D_phidot(a, adot, phi, phidot, m2)
        phi = D_phi(phi, fdot_old)

        a = D_a(a, adot)

        adot = Friedman(a, phi, phidot, m2)

        t -= dt

        tlist.append(t)
        alist.append(a)
    #    adotl.append(adot)
    #    phil.append(phi)
    #    phidotl.append(phidot)

        k += 1

    tar = np.array(tlist)
    arr = np.array(alist)

    return tar, arr


# Make the plots

ax = subplot(111)
subplots_adjust(left=0.25, bottom=0.25)

time, scale = run_dat_shiet(a0, adot0, phi0, phidot0, mass2)

l, = plot(time, scale)

axcolor = 'lightgoldenrodyellow'
axmass = axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
axphi  = axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)

smass = Slider(axmass, 'scalar mass', 0.0, 0.0001, valinit=np.sqrt(mass2))
sfin_phi = Slider(axphi, 'current phi', 0., 1.0, valinit=phi0)

def update(val):
	mass = smass.val
	fin_phi = sfin_phi.val
        newt, news = run_dat_shiet(a0, adot0, fin_phi,
                                   init_phidot(fin_phi, mass**2), mass**2)
        l.set_xdata(newt)
	l.set_ydata(news)
	draw()
smass.on_changed(update)
sfin_phi.on_changed(update)

resetax = axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
def reset(event):
	smass.reset()
	sfin_phi.reset()
button.on_clicked(reset)

show()

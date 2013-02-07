import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def RK4(dt, yn, f, *args):

    k1 = dt * f(*args)
    k2 = dt * f(args[0] + k1/2., *args[1:])
    k3 = dt * f(args[0] + k2/2., *args[1:])
    k4 = dt * f(args[0] + k3, *args[1:])

    return yn + k1/6. + k2/3. + k3/3. + k4/6.


def adv_a(a, phi, m2, H0):

    return np.sqrt((0.3*H0)/a + ((a**2)/6.)*(m2*phi**2))


def adv_phi(phi, a, m2, H0):

    adot = adv_a(a, phi, m2, H0)
    return (1/3.) * (a/adot) * m2 * phi


# initial values

dt = -0.00001    # stepsize
a = 1.
phi = 1.
m2 = 10**(-15)
H0 = np.sqrt((1./6.) * (m2/0.7))

alist = [a]
phil = [phi]
t = 0.
tlist = [0.]
k = 1

while a > 0.0001 and k < 100000:

    phi = RK4(dt, phi, adv_phi, phi, a, m2, H0)

    a = RK4(dt, a, adv_a, a, phi, m2, H0)

    t += dt

    tlist.append(t)
    alist.append(a)
#    adotl.append(adot)
#    phil.append(phi)
#    phidotl.append(phidot)

    k += 1


# Make the plots

tar = np.array(tlist)
arr = np.array(alist)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(tar, arr)

fig.savefig("/home/kilian/public_html/quentin/approx_RK4.png")

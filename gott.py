import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')
import numpy as np
from scipy.integrate import odeint as oink
from scipy.integrate import trapz


def RK4(dt, yn, f, *args):

    k1 = dt * f(*args)
    k2 = dt * f(args[0] + k1/2., *args[1:])
    k3 = dt * f(args[0] + k2/2., *args[1:])
    k4 = dt * f(args[0] + k3, *args[1:])

    return yn + k1/6. + k2/3. + k3/3. + k4/6.


def H_tau(a, exp_int):

    return np.sqrt(0.0001*a**(-2) + 0.2699*a**(-1) + (0.73*a**2)*np.exp(3.*exp_int))


def adv_z(z, H_t, m_star2):

    return - z**2 - 3*H_t*z - m_star2

if __name__=="__main__":

    m_star2 = 1.

    #---------------
    # Gettin Started!
    #-----------------

    # First find a dumb H(tau)

    a0 = 1.0e-20    # make starting size of universe tiny
    t1 = np.arange(0.0, a0*1000, a0)
    t2 = np.arange(t1[-1]+a0, (t1[-1]+a0)*1000, a0*1000) # approx first time range
    t3 = np.arange(t2[-1]+a0*1000, (t2[-1]+a0)*1000, a0*10.0e6) # approx first time range
    t4 = np.arange(t3[-1]+a0*10.0e6, 1.1, 0.001)
    t = np.concatenate([t1, t2, t3, t4])

#    a = oink(first_int, a0, t)
#    a = np.reshape(a, len(a))
    a = np.array([a0])
    for i in range(len(t)-1):
        a = np.append(a, RK4(t[i+1] - t[i], a[i],
                             H_tau, a[i], 0.0) )


    end = np.abs(a - 1.0).argmin()  # index of tau where a = 1
    tau = t[:end]
    tau_0 = t[end]

    print "tau_0 is ", tau_0

    a_dot = H_tau(a, 0.0)

    a = a[:end]         #\
    a_dot = a_dot[:end] # |- first guess stuff
    H = a_dot/a         #/

    print "initial H is: ", H


    # start the plotting

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(tau, a, label=r"$\Lambda$CDM")

    # now to iterate

    difference = 10.  # just to get the while loop rolling
    delta_w0 = 1.0
    old_w = -1.0*np.ones(len(t[:end])) # w for first guess

    k = 0 # let's see how long until convergence

    while (difference > (0.1*delta_w0)):

        zeta = np.array([0.])

        for i in range(len(tau) - 1):
            '''advance all the zeta using old H'''
            zeta = np.append(zeta,
                             RK4(tau[i+1] - tau[i], zeta[i],
                                 adv_z, zeta[i], H[i], m_star2) )

        print "zeta is ", zeta

        w_tau = np.array([])

        for z in zeta:
            '''compute new eqn of state'''
            w_tau = np.append(w_tau, (z**2 - m_star2)/(z**2 + m_star2))

        print "w_tau is ", w_tau

        exp_int = np.array([])

        for i in range(len(H)):
            '''get an array for the exponent in the new H'''
            exp_int = np.append(exp_int, trapz(H[i:]*(1. + w_tau[i:])))

        new_a = np.array([a0])

        for i in range(len(tau)-1):
            new_a = np.append(new_a,
                              RK4(tau[i+1] - tau[i], new_a[i],
                                  H_tau, new_a[i], exp_int[i]) )

        print "new a is", new_a

        new_a_dot = H_tau(new_a, exp_int)

        difference = np.abs(trapz(w_tau) - trapz(old_w))
        delta_w0 = np.abs(1.0 + w_tau[-1])

        old_w = w_tau
        a = new_a
        H = new_a_dot/a

        k += 1


    # last of the plotting

    ax1.plot(tau, a, label=r"QCDM")
    ax1.set_xlabel(r"$\tau$")
    ax1.set_ylabel(r"a")

    ax2 = fig.add_subplot(212)
    ax2.plot(tau, np.abs(1.0-w_tau))
    ax2.set_xlabel(r"$\tau$")
    ax2.set_ylabel(r"$\delta$w($\tau$)")

    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels)

    fig.savefig("/home/kilian/public_html/quentin/gott.png")

    print "the number of runs to converge was k=", k

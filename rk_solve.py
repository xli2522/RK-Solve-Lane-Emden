import numpy as np
import matplotlib.pyplot as plt


def odeModel(state, z, n):
    '''Input: 
        state    initial state vector
        z       independent variable
        n       parameter
        Output:
        deriv Derivatives vector
        '''
    # set initial v and u for each integration step
    v = state[0]
    u = state[1]
    # parameter n=n

    # compute the right side of the ODEs
    dudz = v
    dvdz = -u**n - 2*v/z

    # output
    derive = np.array([dudz, dvdz])
 
    return derive

def rk2(x, t, tau, derivsRK, n):
    '''
    #  Runge-Kutta integrator (2nd order)
    # Input arguments -
    #   x = current value of dependent variable
    #   t = independent variable (usually time)
    #   tau = step size (usually timestep)
    #   derivsRK = right hand side of the ODE; derivsRK is the
    #             name of the function which returns dx/dt
    #             Calling format derivsRK (x,t,param).
    #   param = extra parameters passed to derivsRK
    # Output arguments -
    #   xout = new value of x after a step of size tau
    '''
    half_tau = 0.5*tau
    F1 = derivsRK(x,t, n)  
    t_half = t + half_tau
    xtemp = x + half_tau*F1
    F2 = derivsRK(xtemp,t_half,n) 
    xout = x + tau/2.*(F1 + F2)

    return xout

def main(v0, u0, tau, n):
    '''Returns the righ side of the coupled first-order ODEs; used by RK routines
    Inputs: 
            v0      initial v
            u0      initial u
            tau     integration step size
    Output:
            u, v Results after the integration '''
    u = u0   #initial u
    v = v0   #initial v
    z = 0    #initial z
    uplot = []
    vplot = []
    zplot = []
    
    state = np.array([u, v])
    
    while u >= 0:
        uplot.append(u)
        vplot.append(v)
        zplot.append(z)

        z += tau    #next time step
        # update new state
        state = rk2(state, z, tau, odeModel, n)

        # get updated u, v
        u = state[0]
        v = state[1]

    return  uplot, vplot, zplot

def plot(x, y):
    fig, ax = plt.subplots()
    ax.plot(x, y)
    plt.show()

    return 

u, v, z = main(1, 0.01, 0.01, 2.5)

print(u)
# Solve Lane-Emden with RK2
import numpy as np
import matplotlib.pyplot as plt

'''ODE Model:
    dudxi = v
    dudxi = -u**n - 2*v/xi
    power series:
    dudxi = -u**n - (2*v)+2*v*(xi-1)
    Problem: xi approaches infinity when u approaches 0 for some n
            --> Stopping condition can not be u > 0 alone
    Solution: apply power series'''

tau = 0.01      # appropriate time step
fig, ax = plt.subplots()        # set up plot

# loop through different n values
for n in [1.5, 3, 3.25]:
#for n in range(6):
    xi = 0       # initial xi
    u = 1           # intial u
    v = 0           # initial v
    # container for final solutions of each variable 
    uplot = []
    vplot = []
    xiplot = []
    count = 0       # counter the total number of steps of integration

    # Set stopping conditions: reaches the surface
    while u >= 0:       # only works for small n, ie: n < 4
        # RK2 integration
        half_tau = 0.5*tau
        xi += tau
        if count == 0:
            # use power series approximation when xi = 0 (first two terms)
            v += tau*(-(u+u*half_tau)**n - (2*(v+v*half_tau))+2*(v+v*half_tau)*((xi+xi*half_tau)-1))
        else:
            v += tau*(-(u+u*half_tau)**n - 2*(v+v*half_tau)/(xi+xi*half_tau))
        u += tau*(v+v*half_tau)
        # record solutions
        uplot.append(u) 
        vplot.append(v)
        xiplot.append(xi)
        count += 1
    ax.plot(xiplot, uplot)      # plot xi vs. u under each n

# analytical
for n in [0, 1, 5]:
    uAnalytical = []
    xiAnalytical = []
    xiAnalytical.append()



    
ax.set_xlim([0, 10])        # set x axis limit
ax.set_ylim([0, 1.0])       # set y axis limit
ax.legend(['n = 1.5','n = 3','n = 3.5'])
plt.show()      # show plot

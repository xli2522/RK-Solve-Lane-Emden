# Solve Lane-Emden with RK2
import numpy as np
import matplotlib.pyplot as plt

'''ODE Model:
    dudxi = v
    dudxi = -u**n - 2*v/xi
    Problem: xi approaches infinity when u approaches 0 for some n
            --> Stopping condition can not be u > 0 alone
    Solution: apply power series'''

tau = 0.01      # appropriate time step
fig, ax = plt.subplots()        # set up plot

# loop through different n values
#for n in [1.5, 3, 3.25]:
for n in range(6):
    xi = 0.0001       # initial xi
    u = 1           # intial u
    v = 0           # initial v
    # container for final solutions of each variable 
    uplot = []
    vplot = []
    xiplot = []
    count = 0       # counter the total number of steps of integration

    # Set stopping conditions: reaches the surface
    while u >= 0 and xi < 30:
        # add initial conditions
        half_tau = 0.5*tau
        xi += tau
        
        # RK2 integration
        v += tau*(-(u+u*half_tau)**n - 2*(v+v*half_tau)/(xi+xi*half_tau))
        u += tau*(v+v*half_tau)
        # record solutions
        uplot.append(u) 
        vplot.append(v)
        xiplot.append(xi)
        count += 1
    ax.plot(xiplot, uplot)      # plot xi vs. u under each n

ax.set_xlim([0, 10])        # set x axis limit
ax.set_ylim([0, 1.0])       # set y axis limit
plt.show()      # show plot

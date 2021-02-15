import numpy as np
import matplotlib.pyplot as plt

tau = 0.001      # time step 
fig, ax = plt.subplots()
# try different n values
for n in [1.5, 3, 3.25]:
    xi = 0.0001       # initial condition for xi
    u = 1           # intial condition for u, theta, the nondimensionalized radius
    v = 0
    # collectors for our final solutions for each variable 
    uplot = []
    vplot = []
    xiplot = []
    count = 0
    while u >= 0 and xi < 20:
        # add initial conditions
        half_tau = 0.5*tau
        xi += half_tau
        v += tau*(-(u+u*half_tau)**n - 2*(v+v*half_tau)/(xi+xi*half_tau))
        u += tau*(v+v*half_tau)

        uplot.append(u)
        vplot.append(v)
        xiplot.append(xi)
        count += 1
    ax.plot(xiplot, uplot)
plt.show()


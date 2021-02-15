# sample solution using scipy's odeint
# reference stackoverflow

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

def poly(y,xi,n):  
    theta, phi = y  
    dydxi = [-phi/(xi**2), (xi**2)*theta**n]  
    return dydxi

fig, ax = plt.subplots()

y0 = [1.,0.]  
xi = np.linspace(10e-4, 16., 201)

for n in range(6):
    sol = odeint(poly, y0, xi, args=(n,))
    ax.plot(xi, sol[:, 0])

ax.axhline(y = 0.0, linestyle = '--', color = 'k')
ax.set_xlim([0, 10])
ax.set_ylim([0, 1.0])
ax.set_xlabel(r"$\xi$")
ax.set_ylabel(r"$\theta$")
plt.show()
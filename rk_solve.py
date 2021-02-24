# Solve Lane-Emden with RK2
import numpy as np
import matplotlib.pyplot as plt
import math

'''ODE Model:
    dudxi = v
    dudxi = -u**n - 2*v/xi
    power series:
    dudxi = -u**n - (2*v)+2*v*(xi-1)
    Problem: xi approaches infinity when u approaches 0 for some n
            --> Stopping condition can not be u > 0 alone
    Solution: apply power series'''

tau = 0.01      # appropriate time step
rho_c=125.09001252015766 # this is the calculated central density using our model, calculated below but hard-coded here for plotting purposes for Question 3
fig, ax = plt.subplots()        # set up plot for Question 1

fig2, ax2 =  plt.subplots() 	# set up the plot for Question 3
thetas3point25=[] 	# this is a list that will hold the values of theta^3.25 for plotting
zetas3point25=[] 	# this is a list that will hold the values of zeta for plotting (to calculate r/R_sun)

#hard-code in the values of rho and r from the solar model on OWL
r_over_R_sun=[0.007,0.02,0.09,0.22,0.32,0.42,0.52,0.60,0.71,0.81,0.91,0.96,0.99,0.995,0.999,1.000]
rho=[150,146,95.73,28.72,9.77,3.22,1.05,0.500,0.177,0.0766,0.0194,4.85e-3,2.56e-4,4.83e-5,1.29e-6,2.18e-7]

ax2.plot(r_over_R_sun, rho, label="Standard Solar Model") # plot the standard solar model values

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
    rho_ours=[] # this is just for Question 3
    zetas3point25=[]
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

        if n==3.25:
            rho_ours.append(rho_c*(u**3.25)) # get the appropriate value for t
            zetas3point25.append(xi)

    # if n==3.25:
    # 	zetas3point25=xiplot
    # 	thetas3point25=uplot


    ax.plot(xiplot, uplot, label="n="+str(n))      # plot xi vs. u under each n
    
#define functions for the analytical solutions for n=0,1,5
def theta0(zeta):
	theta=1-(zeta**2/6)
	return theta
def theta1(zeta):
	theta=math.sin(zeta)/zeta
	return theta
def theta5(zeta):
	theta=(1+(zeta**2/3))**(-0.5)
	return theta

#generate the analytical solutions to be plotted
thetas_0=[]
thetas_1=[]
thetas_5=[]
for i in range(len(xiplot)):
	thetas_0.append(theta0(xiplot[i]))
	thetas_1.append(theta1(xiplot[i]))
	thetas_5.append(theta5(xiplot[i]))

#plot the analytic solutions over top of the numeric solutions
ax.plot(xiplot,thetas_0, label="n=0 (analytic)")
ax.plot(xiplot,thetas_1, label="n=1 (analytic)")
ax.plot(xiplot,thetas_5, label="n=5 (analytic)")

# Calculation of q2 (a)
xiRe = xiplot.pop()      # xi at the surface        7.979999999999874
vRe = vplot.pop()       # dudxi at the surface      -0.02996372129871972
Rsun = 6.96*10**10      # cm
Msun = 1.99*10**33      # g
rhoc = - (Msun * xiRe)/(4*np.pi*Rsun**3*vRe)
print(rhoc)         # 125.09001252015766

ax.set_xlim([0, max(xiplot)])        # set x axis limit
ax.set_ylim([0, 1.0])       # set y axis limit
ax.set_xlabel('\u03B6')     # set axis titles
ax.set_ylabel('\u03B8')
ax.set_title("Numeric and Analytic Solutions to the Lane-Emden Equation")       # set plot title
ax.legend()

#now do the plotting for Question 3
# rho_ours=[] # initialize a list that will hold the values of rho
# for i in range(len(thetas3point25)):
# 	rho.append(rhoc*thetas3point25[i]**3.25)
# # print(len(rho),len(zetas3point25))
ax2.plot(zetas3point25, rho_ours, label="from our 3.25 solution")
ax2.legend()
ax2.set_title(" Density "+'\u03C1'+ " as a Function of Radius")
ax2.set_xlabel("Radius [solar radii]")
ax2.set_ylabel("Density [g cm$^{-3}$]")



plt.show()      # show plot

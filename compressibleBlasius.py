#=====================================================
# Department of Mechanical Engineering, Purdue University
# ME 614: Computational Fluid Dynamics
# Fall 2018
# Compressible Blasius Boundary Layer Solver
# Julien Brillon
# Python 2.7.15
#=====================================================
# Import libraries
import numpy as np # NumPy: contains basic numerical routines
import scipy # SciPy: contains additional numerical routines to numpy
import matplotlib.pyplot as plt # Matlab-like plotting
import scipy.sparse as scysparse
import scipy.sparse.linalg
from scipy.optimize import fsolve, bisect, minimize
from scipy import interpolate
#=====================================================
data_fileType = 'txt'
subdirectories = ['Data/','Figures/']
#=====================================================
def sysFG(eta,u):
    global Pr, gamma, Me, omega, theta, viscosity_model
    # Initial conditions
    f = np.array([u[0],u[1],u[2]]) # f, f', f''
    g = np.array([u[3],u[4]]) # g, g'
    
    T_ratio = g[0] + 0.5*(Me**2.0)*(gamma-1.0)*(g[0]-f[1]**2.0)
    
    if viscosity_model == 'Sutherland':
        invC = (T_ratio**(-0.5))*((T_ratio+theta)/(1.0+theta))
    elif viscosity_model == 'Power law':
        invC = T_ratio**(1.0-omega) # Power law

    K = ((2.0*(gamma-1.0)*Me**2.0)/(2.0+(gamma-1.0)*Me**2.0)) # coefficient ue^2/(h0)e^2 in the energy equation

    # System of 1st order ODEs
    df1 = f[1]
    df2 = f[2]
    df3 = -invC*f[0]*f[2]
    dg1 = g[1]
    dg2 = -Pr*invC*f[0]*g[1] - K*(Pr-1.0)*f[2]*(f[2]-invC*f[0]*f[1])

    return np.array([df1,df2,df3,dg1,dg2])
#=====================================================
def RK4_FG(y):
    global n, h, eta

    for i in range(0,n):
        k1=h*sysFG(eta[i], y[i,:])
        k2=h*sysFG(eta[i]+0.5*h, y[i,:]+0.5*k1)
        k3=h*sysFG(eta[i]+0.5*h, y[i,:]+0.5*k2)
        k4=h*sysFG(eta[i]+h, y[i,:]+k3)
        y[i+1,:]=y[i,:] + (1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)

    return y
#=====================================================
def nonlinear_FG(solverInput):
    global n, m, f00, f10, f1inf, g0inf, wall_temperature_BC_type

    if wall_temperature_BC_type == 'Dirichlet':
        global g00
        f20, g10 = solverInput
    elif wall_temperature_BC_type == 'Neumann':
        global g10
        f20, g00 = solverInput
    
    u = np.zeros((n+1,m+2))
    u[0,:] = np.array([f00,f10,f20,g00,g10])

    u = RK4_FG(u)

    if wall_temperature_BC_type == 'Dirichlet':
        return (u[n,1] - f1inf, u[n,3] - g0inf)
    elif wall_temperature_BC_type == 'Neumann':
        return (u[n,1] - f1inf, u[n,3] - g0inf)
#=====================================================
def update_FG(f20,g0):
    global n, m, f00, f10, wall_temperature_BC_type
    
    if wall_temperature_BC_type == 'Dirichlet':
        global g00
        g10 = g0
    elif wall_temperature_BC_type == 'Neumann':
        global g10
        g00 = g0

    u = np.zeros((n+1,m+2))
    u[0,:] = np.array([f00,f10,f20,g00,g10])
    u = RK4_FG(u)

    f = np.empty((n+1,m)) # f, f', f''
    for i in range(0,m):
        f[:,i] = u[:,i]
    
    g = np.empty((n+1,m-1)) # g, g'
    for i in range(m,m+2):
        g[:,i-m] = u[:,i]

    return (f, g)
#=====================================================
def plot_F(f):
    global eta
    # Post processing
    figure_title = "Velocity Profile"
    figure_title_print = figure_title
    print('Plotting: ' + figure_title_print)
    fig, ax = plt.subplots(figsize=(7,7))
    ax.set_title(figure_title)
    ax.set_ylabel(r"$\frac{y}{x}\sqrt{Re_{x}}$")
    ax.set_xlabel(r"$u/u_{e}$")
    
    ax.set_xlim(0.0, 1.5)
    # ax.set_ylim(0.0, 40)
    # ax.set_aspect(0.06)
    # plt.ylim(a,b)
    plt.plot(f[:,0],eta,label="f")
    plt.plot(f[:,1],eta,label="f'")
    plt.plot(f[:,2],eta,label="f''")

    plt.grid()
    plt.tight_layout()
    leg = plt.legend(loc='best', ncol=1, shadow=True, fancybox=True, fontsize=8)
    
    plt.show()
#=====================================================
def plot_G(g):
    global eta
    # Post processing
    figure_title = "Total Enthalpy Profile"
    figure_title_print = figure_title
    print('Plotting: ' + figure_title_print)
    fig, ax = plt.subplots()
    ax.set_title(figure_title)
    ax.set_ylabel(r"$\frac{y}{x}\sqrt{Re_{x}}$")
    ax.set_xlabel(r"$h_0/(h_0)_{e}$")

    # plt.ylim(a,b)

    plt.plot(g[:,0],eta,label="g")
    plt.plot(g[:,1],eta,label="g'")

    plt.grid()
    plt.tight_layout()
    leg = plt.legend(loc='best', ncol=1, shadow=True, fancybox=True, fontsize=8)

    plt.show()
#=====================================================
def plot_T(T_ratio):
    global eta, y_var
    # Post processing
    figure_title = "Temperature Profile"
    figure_title_print = figure_title
    print('Plotting: ' + figure_title_print)
    fig, ax = plt.subplots(figsize=(7,7))
    ax.set_title(figure_title)

    if y_var == 'eta':
        ax.set_ylabel(r"$\eta$")
    elif y_var == 'y_xRe_x':
        ax.set_ylabel(r"$\frac{y}{x}\sqrt{Re_{x}}$")
    ax.set_xlabel(r"$T/T_{e}$")
    
    # ax.set_xlim(0.0, 4.0) # 4.0
    # ax.set_ylim(0.0, 10) #40.0
    # ax.set_aspect(0.5)

    plt.plot(T_ratio,eta)

    plt.grid()
    plt.tight_layout()

    plt.show()
#=====================================================
global Pr, n, m, h, eta, f00, f10, f1inf, g0inf, gamma, Me, omega, theta
global viscosity_model, wall_temperature_BC_type, y_var
viscosity_models = ['Sutherland','Power law']
wall_temperature_BC_types = ['Dirichlet','Neumann']
a = 0.0
b = 6.0
h = 0.01
n = (b-a)/h
n = np.int(n)
m = 3
eta = np.linspace(a,b,n+1)
# Flow parameters
Pr = 0.75 # Prandtl number
gamma = 1.4
omega = 0.76 # Power law parameter for viscosity [Driest 1952]
theta = 0.505 # Sutherland non-dimensional S/Te [Driest 1952]
viscosity_model = viscosity_models[0]
wall_temperature_BC_type = wall_temperature_BC_types[0]
# Boundary conditions
f00 = 0.0
f10 = 0.0
f1inf = 1.0
g0inf = 1.0

# Specify mach number
Me = raw_input('Specify mach number: ')
Me = np.float(Me)
# Me = 10.0

subdirectory = subdirectories[0] + wall_temperature_BC_type + '/'
guess_store_name = 'smart_guess_list_' + wall_temperature_BC_type +'.'+ data_fileType
Me_study, f20_guess_store, g0_guess_store = np.loadtxt(subdirectory+guess_store_name,unpack=True)

if wall_temperature_BC_type == 'Dirichlet':
    global g00
    T_ratio_wall = 1.0 # Tw/Te
    g00 = (T_ratio_wall+0.5*(Me**2.0)*(gamma-1.0)*(f10**2.0))/(1.0+0.5*(Me**2.0)*(gamma-1.0))
    # Me_study = np.array([0.0,2.0,3.0,4.0,6.0,7.0,8.0,8.5,9.0,9.25,9.5,9.75,9.875,10.0])
    # f20_guess_store = np.array([0.4696,0.4746,0.4810,0.4900,0.5145,0.5291,0.5449,0.5530,0.5613,0.5655,0.5698,0.5740,0.5761,0.5783])
    # g0_guess_store = np.array([0.0,0.1765,0.2589,0.3130,0.3800,0.4046,0.4266,0.4370,0.4470,0.4519,0.4568,0.4616,0.4640,0.4663])

elif wall_temperature_BC_type == 'Neumann':
    global g10
    g10 = 0.0 # adiabatic wall
    # Me_study = np.array([0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,17.0,18.0,19.0,20.0])
    # f20_guess_store = np.array([0.4696,0.4844,0.5263,0.5818,0.6408,0.6990,0.7551,0.8086,0.8597,0.8843,0.9084,0.9320,0.9550])
    # g0_guess_store = np.array([1.0000,0.9692,0.9572,0.9623,0.9699,0.9766,0.9818,0.9858,0.9889,0.9902,0.9914,0.9924,0.9933])

# Initial guesses (magical for adiabatic (Neumann))
if Me > np.max(Me_study):
    # Extrapolate guess for f20
    f20_guess_int = np.poly1d(np.polyfit(Me_study,f20_guess_store,2))
    f20_guess = f20_guess_int(Me)
    # Extrapolate guess for g00
    g0_guess_int = np.poly1d(np.polyfit(Me_study,g0_guess_store,2))
    g0_guess = g0_guess_int(Me)
    # Open file to store the correct guess after convergence
    subdirectory = subdirectories[0] + wall_temperature_BC_type + '/'
    guess_file = open(subdirectory + guess_store_name,"a+")
    # guess_file.write("%.2f %.4f %.4f" % (25.0,7.0,1.0))
    # guess_file.close()
else:
    mach_index = list(Me_study).index(Me)
    f20_guess = f20_guess_store[mach_index]
    g0_guess = g0_guess_store[mach_index]

# Exactly solve the coupled nonlinear ODEs
if wall_temperature_BC_type == 'Dirichlet':
    # for g0_guess in np.arange(0.5,3.5,0.5):
    f20, g10 = fsolve(nonlinear_FG, (f20_guess,g0_guess),xtol=1.0e-10)
    f, g = update_FG(f20,g10)
    if Me > np.max(Me_study):
        guess_file.write("%.4f %.4f %.4f\r\n" % (Me,f20,g10))
        guess_file.close()
    # if (np.abs(g[n,0] - g0inf) < 1.0e-10) and (np.abs(f[n,1] - f1inf) < 1.0e-10):
    #   print('END')
    #   break
elif wall_temperature_BC_type == 'Neumann':
    f20, g00 = fsolve(nonlinear_FG, (f20_guess,g0_guess),xtol=1.0e-10)
    f, g = update_FG(f20,g00)
    if Me > np.max(Me_study):
        guess_file.write("%.4f %.4f %.4f\r\n" % (Me,f20,g00))
        guess_file.close()

# Print error to check that condition is met
print("f20_error = %.3e" % (f[n,1]-f1inf))
print("g00_error = %.3e" %(g[n,0]-g0inf))

# Post processing
T_ratio = g[:,0] + 0.5*(Me**2.0)*(gamma-1.0)*(g[:,0]-(f[:,1])**2.0)

rho_ratio = np.ones(n+1)
for i in range(0,n+1):
    rho_ratio[i] = np.trapz(T_ratio[0:i+1])

y_xRe_x = np.sqrt(2.0)*rho_ratio

y_var = 'eta'

print('Edge mach number, Me = %.4f' % Me)
# plot_F(f)
# plot_G(g)
# plot_T(T_ratio)

# ---------------------------------------------------------------
#                       Export data 
# ---------------------------------------------------------------
Me_select = np.loadtxt('Me_select.txt')
Me_select = Me_select.tolist()
if Me in Me_select:
    subdirectory = subdirectories[0] + wall_temperature_BC_type + '/'
    filename = "eta_M_%i_BC_%s" % (np.int(Me),wall_temperature_BC_type)
    np.savetxt(subdirectory+filename+'.'+data_fileType, eta)
    filename = "f_M_%i_BC_%s" % (np.int(Me),wall_temperature_BC_type)
    np.savetxt(subdirectory+filename+'.'+data_fileType, f)
    filename = "g_M_%i_BC_%s" % (np.int(Me),wall_temperature_BC_type)
    np.savetxt(subdirectory+filename+'.'+data_fileType, g)
    filename = "T_ratio_M_%i_BC_%s" % (np.int(Me),wall_temperature_BC_type)
    np.savetxt(subdirectory+filename+'.'+data_fileType, T_ratio)
    filename = "yxRex_M_%i_BC_%s" % (np.int(Me),wall_temperature_BC_type)
    np.savetxt(subdirectory+filename+'.'+data_fileType, y_xRe_x)

# y_var = 'y_xRe_x'
# eta = y_xRe_x # for plotting
# plot_F(f)
# plot_G(g)
# plot_T(T_ratio)

print('-----------------------------------------------------')
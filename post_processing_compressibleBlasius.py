#=====================================================
# Department of Mechanical Engineering, Purdue University
# ME 614: Computational Fluid Dynamics
# Fall 2018
# Julien Brillon
# Python 2.7.15
#=====================================================
# Import libraries
import numpy as np # NumPy: contains basic numerical routines
import scipy # SciPy: contains additional numerical routines to numpy
import matplotlib.pyplot as plt # Matlab-like plotting
import scipy.sparse as scysparse
import scipy.sparse.linalg
from scipy.optimize import fsolve
#=====================================================
data_fileType = 'txt'
subdirectories = ['Data/','Figures/']
#=====================================================
# Me_study = [0,2,4,6,8,10,12,16,20]
Me_select = np.loadtxt('Me_select.txt')
Me_select = Me_select.tolist()
# Me_study = [0,2,4,8,10]
wall_temperature_BC_types = ['Dirichlet','Neumann']
wall_temperature_BC_type = wall_temperature_BC_types[0]

if wall_temperature_BC_type == 'Dirichlet':
    figure_title = "Velocity Profiles in a Laminar, Compressible\n Boundary Layer over an Isothermal Flat Plate"
    figure_title_print = "Velocity Profiles in a Laminar, Compressible Boundary Layer over an Isothermal Flat Plate"
elif wall_temperature_BC_type == 'Neumann':
    figure_title = "Velocity Profiles in a Laminar, Compressible\n Boundary Layer over an Adiabatic Flat Plate"
    figure_title_print = "Velocity Profiles in a Laminar, Compressible Boundary Layer over an Adiabatic Flat Plate"

print('Plotting: ' + figure_title_print)
fig, ax = plt.subplots(figsize=(7,7))
ax.set_title(figure_title)
ax.set_ylabel(r"$\eta$")
ax.set_xlabel(r"$u/u_{e}$")
# ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0,6.0)
ax.set_aspect(0.25)

for Me in Me_select:
    subdirectory=subdirectories[0] + wall_temperature_BC_type + '/'
    filename = "eta_M_%i_BC_%s" % (np.int(Me),wall_temperature_BC_type)
    eta = np.loadtxt(subdirectory+filename+'.'+data_fileType,unpack=True)
    filename = "f_M_%i_BC_%s" % (np.int(Me),wall_temperature_BC_type)
    f = np.loadtxt(subdirectory+filename+'.'+data_fileType,unpack=False)
    name = r'$M_{e} = %i$' % np.int(Me)
    plt.plot(f[:,1],eta,label=name)

plt.grid()
plt.tight_layout()
leg = plt.legend(loc='best', ncol=1, shadow=True, fancybox=True, fontsize=8)
print(' ... Saving figure ...')
figure_name = "velocity_profile_%s" % wall_temperature_BC_type
figure_name = figure_name
figure_fileType = 'png'
subdirectory = subdirectories[1]
plt.savefig(subdirectory + figure_name + '.' + figure_fileType,format=figure_fileType,dpi=500)
plt.close()

if wall_temperature_BC_type == 'Dirichlet':
    figure_title = "Temperature Profiles in a Laminar, Compressible\n Boundary Layer over an Isothermal Flat Plate"
    figure_title_print = "Temperature Profiles in a Laminar, Compressible Boundary Layer over an Isothermal Flat Plate"
elif wall_temperature_BC_type == 'Neumann':
    figure_title = "Temperature Profiles in a Laminar, Compressible\n Boundary Layer over an Adiabatic Flat Plate"
    figure_title_print = "Temperature Profiles in a Laminar, Compressible Boundary Layer over an Adiabatic Flat Plate"

print('Plotting: ' + figure_title_print)
fig, ax = plt.subplots(figsize=(7,7))
ax.set_title(figure_title)
ax.set_ylabel(r"$\eta$")
ax.set_xlabel(r"$T/T_{e}$")
# ax.set_xlim(0.0,.0)
ax.set_ylim(0.0,6.0)

if wall_temperature_BC_type == 'Dirichlet':
    ax.set_aspect(1.5)
elif wall_temperature_BC_type == 'Neumann':
    ax.set_aspect(15.0)

for Me in Me_select:
    subdirectory=subdirectories[0] + wall_temperature_BC_type + '/'
    filename = "eta_M_%i_BC_%s" % (np.int(Me),wall_temperature_BC_type)
    eta = np.loadtxt(subdirectory+filename+'.'+data_fileType,unpack=True)
    filename = "T_ratio_M_%i_BC_%s" % (np.int(Me),wall_temperature_BC_type)
    T_ratio = np.loadtxt(subdirectory+filename+'.'+data_fileType,unpack=True)
    name = r'$M_{e} = %i$' % np.int(Me)
    plt.plot(T_ratio,eta,label=name)

plt.grid()
plt.tight_layout()
leg = plt.legend(loc='best', ncol=1, shadow=True, fancybox=True, fontsize=8)
print(' ... Saving figure ...')
figure_name = "temperature_profile_%s" % wall_temperature_BC_type
figure_name = figure_name
figure_fileType = 'png'
subdirectory = subdirectories[1]
plt.savefig(subdirectory + figure_name + '.' + figure_fileType,format=figure_fileType,dpi=500)
plt.close()

# filename = "g_M_%i_BC_%s" % (np.int(Me),wall_temperature_BC_type)
# np.savetxt(subdirectory+filename+'.'+data_fileType, g)
# filename = "T_ratio_M_%i_BC_%s" % (np.int(Me),wall_temperature_BC_type)
# np.savetxt(subdirectory+filename+'.'+data_fileType, T_ratio)
# filename = "yxRex_M_%i_BC_%s" % (np.int(Me),wall_temperature_BC_type)
# np.savetxt(subdirectory+filename+'.'+data_fileType, y_xRe_x)



# #=====================================================
# def plot_F():
#   global eta
#   # Post processing
#   figure_title = "Velocity Profile"
#   figure_title_print = figure_title
#   print('Plotting: ' + figure_title_print)
#   fig, ax = plt.subplots(figsize=(7,7))
#   ax.set_title(figure_title)
#   ax.set_ylabel(r"$\frac{y}{x}\sqrt{Re_{x}}$")
#   ax.set_xlabel(r"$u/u_{e}$")
    
#   ax.set_xlim(0.0, 1.0)
#   # ax.set_ylim(0.0, 40)
#   # ax.set_aspect(0.06)
#   # plt.ylim(a,b)
#   plt.plot(f[:,0],eta,label="f")
#   plt.plot(f[:,1],eta,label="f'")
#   plt.plot(f[:,2],eta,label="f''")

#   plt.grid()
#   plt.tight_layout()
#   leg = plt.legend(loc='best', ncol=1, shadow=True, fancybox=True, fontsize=8)
    
#   plt.show()
    
#   # print(' ... Saving figure ...')
#   # figure_name = "iters_vs_omega_n_%i_log10Re_%i" % (n[j],np.abs(Re))
#   # figure_name = figure_name+file_ext
#   # figure_fileType = 'png'
#   # subdirectory = subdirectories[1]
#   # plt.savefig(subdirectory + figure_name + '.' + figure_fileType,format=figure_fileType,dpi=500)
#   # plt.close()
# #=====================================================
# def plot_G(g):
#   global eta
#   # Post processing
#   figure_title = "Total Enthalpy Profile"
#   figure_title_print = figure_title
#   print('Plotting: ' + figure_title_print)
#   fig, ax = plt.subplots()
#   ax.set_title(figure_title)
#   ax.set_ylabel(r"$\frac{y}{x}\sqrt{Re_{x}}$")
#   ax.set_xlabel(r"$h_0/(h_0)_{e}$")

#   # plt.ylim(a,b)

#   plt.plot(g[:,0],eta,label="g")
#   plt.plot(g[:,1],eta,label="g'")

#   plt.grid()
#   plt.tight_layout()
#   leg = plt.legend(loc='best', ncol=1, shadow=True, fancybox=True, fontsize=8)

#   plt.show()
#   # leg = plt.legend(loc='best', ncol=1, shadow=True, fancybox=True, fontsize=8)
#   # print(' ... Saving figure ...')
#   # figure_name = "iters_vs_omega_n_%i_log10Re_%i" % (n[j],np.abs(Re))
#   # figure_name = figure_name+file_ext
#   # figure_fileType = 'png'
#   # subdirectory = subdirectories[1]
#   # plt.savefig(subdirectory + figure_name + '.' + figure_fileType,format=figure_fileType,dpi=500)
#   # plt.close()
# #=====================================================
# def plot_T(T_ratio):
#   global eta, y_var
#   # Post processing
#   figure_title = "Temperature Profile"
#   figure_title_print = figure_title
#   print('Plotting: ' + figure_title_print)
#   fig, ax = plt.subplots(figsize=(7,7))
#   ax.set_title(figure_title)

#   if y_var == 'eta':
#       ax.set_ylabel(r"$\eta$")
#   elif y_var == 'y_xRe_x':
#       ax.set_ylabel(r"$\frac{y}{x}\sqrt{Re_{x}}$")
#   ax.set_xlabel(r"$T/T_{e}$")
    
#   # ax.set_xlim(0.0, 4.0) # 4.0
#   # ax.set_ylim(0.0, 10) #40.0
#   # ax.set_aspect(0.5)

#   plt.plot(T_ratio,eta)

#   plt.grid()
#   plt.tight_layout()

#   plt.show()
#   # leg = plt.legend(loc='best', ncol=1, shadow=True, fancybox=True, fontsize=8)
#   # print(' ... Saving figure ...')
#   # figure_name = "iters_vs_omega_n_%i_log10Re_%i" % (n[j],np.abs(Re))
#   # figure_name = figure_name+file_ext
#   # figure_fileType = 'png'
#   # subdirectory = subdirectories[1]
#   # plt.savefig(subdirectory + figure_name + '.' + figure_fileType,format=figure_fileType,dpi=500)
#   # plt.close()
# #=====================================================
from PyPonding.structures import ElasticBeam2d
from math import pi, cos, cosh
import numpy as np
import matplotlib.pyplot as plt
import os

# Define units
inch    = 1.0
kip     = 1.0
lb      = kip/1000.0
ft      = 12.0*inch
ksi     = kip/inch**2
psf     = lb/ft**2
pcf     = psf/ft
in_per_ft = inch/ft
m       = inch/0.0254
kN      = kip/4.4482216153

# Input parameters
L       = 40*ft         # Beam span
S       = 5*ft          # Tributary width
E       = 29000.0*ksi   # Modulus of elasticity
I       = 75*inch**4    # Moment of inertia
gamma   = 62.4*pcf      # Unit weight of water
qD      = 0*psf         # Dead load
slope   = 0.25*in_per_ft # Beam slope

zw_over_zh = 0.5
zh      = slope*L
zw      = zw_over_zh*zh # Elevation of water


# Build beam object
beam = ElasticBeam2d(L,S,E,I,gamma,zj=slope*L,qD=qD)
beam.num_elements = 60

# Run Ponding Analysis
resultsP = beam.run_analysis_OPS('IterativeLevel',target_zw=zw)
VmaxP = np.amax(np.absolute(resultsP.shear_along_length))
MmaxP = np.amax(resultsP.bending_moment_along_length)
x_at_MmaxP = resultsP.position_along_length[np.argmax(resultsP.bending_moment_along_length)][0]
TotalLoadP = resultsP.shear_along_length[0]-resultsP.shear_along_length[-1]

# beam.plot_deformed(zw=zw)

# Run First-Order Analysis
beam.include_ponding_effect = False
results1 = beam.run_analysis_OPS('IterativeLevel',target_zw=zw)
Vmax1 = np.amax(np.absolute(results1.shear_along_length))
Mmax1 = np.amax(results1.bending_moment_along_length)
x_at_Mmax1 = results1.position_along_length[np.argmax(results1.bending_moment_along_length)][0]
TotalLoad1 = results1.shear_along_length[0]-results1.shear_along_length[-1]


# Print Results
Cs = (gamma*S*L**4)/(pi**4*E*I)
print('\nCs = %.3f' % Cs)
print('Basic Amplification 1/(1-Cs) = %.3f' % (1/(1-Cs)))
x = pi/2*Cs**0.25
print('More Exact Amplification     = %.3f\n' % ((1/cos(x)-1/cosh(x))/x**2)) # from Silver (2010)
print('                      with ponding           no ponding          amplification')
print('Total Load          %9.3f kip         %9.3f kip          %9.3f' % (TotalLoadP,TotalLoad1,TotalLoadP/TotalLoad1))
print('Max Shear           %9.3f kip         %9.3f kip          %9.3f' % (VmaxP,Vmax1,VmaxP/Vmax1))
print('Max Moment          %9.3f kip-in      %9.3f kip-in       %9.3f' % (MmaxP,Mmax1,MmaxP/Mmax1))
print('Pos. of Mmax      %5.1f in (%4.3f*L)    %5.1f in (%4.3f*L)\n' % (x_at_MmaxP,x_at_MmaxP/L,x_at_Mmax1,x_at_Mmax1/L))


# Plot Results
make_basic_plots = False
make_fancy_plots_US = False
make_fancy_plots_Metric = True

save_folder = 'figures'
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

plt.rc('axes',labelsize=8)
plt.rc('axes',titlesize=8)
plt.rc('legend',fontsize=8)
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)
#plt.rcParams["text.usetex"] = False


if make_basic_plots:
    plt.figure()
    plt.plot(resultsP.position_along_length,resultsP.bending_moment_along_length,label='With Ponding')
    plt.plot(results1.position_along_length,results1.bending_moment_along_length,label='No Ponding')
    plt.legend()
    plt.xlabel('Position along length of beam (in.)')
    plt.ylabel('Bending Moment (kip-in.)')
    plt.xlim([0,L])

    plt.figure()
    plt.plot(resultsP.position_along_length,resultsP.shear_along_length,label='With Ponding')
    plt.plot(results1.position_along_length,results1.shear_along_length,label='No Ponding')
    plt.legend()
    plt.xlabel('Position along length of beam (in.)')
    plt.ylabel('Shear (kips)')
    plt.xlim([0,L])

    plt.figure()
    plt.plot(resultsP.position_along_length,resultsP.bending_moment_along_length/results1.bending_moment_along_length,label='Bending Moment')
    plt.plot(results1.position_along_length,resultsP.shear_along_length/results1.shear_along_length,label='Shear')
    plt.legend()
    plt.xlabel('Position along length of beam (in.)')
    plt.ylabel('Amplification Factor')
    plt.xlim([0,L])
    #plt.ylim([0,2*(1/(1-Cs))])


if make_fancy_plots_US:
    # Moment Diagram
    plt.figure(figsize=(3.5,2.5))
    plt.axes(position=[0.17,0.17,0.80,0.80])
    plt.plot(resultsP.position_along_length,resultsP.bending_moment_along_length,label='With Ponding')
    plt.plot(results1.position_along_length,results1.bending_moment_along_length,label='No Ponding',color='tab:grey')
    plt.plot([x_at_MmaxP,x_at_Mmax1],[MmaxP,Mmax1],linestyle='None',marker='.',color='k')
    plt.plot([0,x_at_MmaxP],[MmaxP,MmaxP],'--',color='k')
    plt.plot([0,x_at_Mmax1],[Mmax1,Mmax1],'--',color='k')
    plt.text(4,MmaxP-7,f'{MmaxP:.03f}', fontsize=8)
    plt.text(79,Mmax1+2.5,f'{Mmax1:.03f}', fontsize=8)
    plt.legend()
    plt.xlabel('Position along length of beam (in.)')
    plt.ylabel('Bending Moment (kip-in.)')
    plt.xlim([0,L])
    plt.ylim(0,100)
    plt.text(35,10,r'$AF_{Moment}=\dfrac{%.03f \:\mathrm{kip-in.}}{%.03f \:\mathrm{kip-in.}} = %.03f$' % (MmaxP,Mmax1,MmaxP/Mmax1),fontsize=8)

    # Shear Diagram
    plt.figure(figsize=(3.5,2.5))
    plt.axes(position=[0.17,0.17,0.80,0.80])
    plt.plot(resultsP.position_along_length,resultsP.shear_along_length,label='With Ponding')
    plt.plot(results1.position_along_length,results1.shear_along_length,label='No Ponding',color='tab:grey')
    plt.plot([0,0], [VmaxP,Vmax1,],linestyle='None',marker='.',color='k')
    plt.plot([0,150],[VmaxP,VmaxP],'--',color='k')
    plt.plot([0,150],[Vmax1,Vmax1],'--',color='k')
    plt.text(95,VmaxP-0.13,f'{VmaxP:.03f}', fontsize=8)
    plt.text(95,Vmax1-0.13,f'{Vmax1:.03f}', fontsize=8)
    plt.legend()
    plt.xlabel('Position along length of beam (in.)')
    plt.ylabel('Shear (kips)')
    plt.xlim([0,L])
    plt.ylim([-0.4,1.4])
    plt.yticks(np.arange(-0.4,1.6,0.2))
    plt.text(200,0.5,r'$AF_{Shear}=\dfrac{%.03f \:\mathrm{kips}}{%.03f \:\mathrm{kips}} = %.03f$' % (VmaxP,Vmax1,VmaxP/Vmax1),fontsize=8)


if make_fancy_plots_Metric:
    # Moment Diagram
    plt.figure(figsize=(3.5,2.5))
    plt.axes(position=[0.15,0.17,0.82,0.80])
    plt.plot(resultsP.position_along_length/m,resultsP.bending_moment_along_length/(kN*m),label='With Ponding')
    plt.plot(results1.position_along_length/m,results1.bending_moment_along_length/(kN*m),label='No Ponding',color='tab:grey')
    plt.plot([x_at_MmaxP/m,x_at_Mmax1/m],[MmaxP/(kN*m),Mmax1/(kN*m)],linestyle='None',marker='.',color='k')
    plt.plot([0,x_at_MmaxP/m],[MmaxP/(kN*m),MmaxP/(kN*m)],'--',color='k')
    plt.plot([0,x_at_Mmax1/m],[Mmax1/(kN*m),Mmax1/(kN*m)],'--',color='k')
    plt.text(4/m,(MmaxP-7)/(kN*m),f'{MmaxP/(kN*m):.03f}', fontsize=8)
    plt.text(79/m,(Mmax1+2.5)/(kN*m),f'{Mmax1/(kN*m):.03f}', fontsize=8)
    plt.legend()
    plt.xlabel('Position along length of beam (m)')
    plt.ylabel('Bending Moment (kN-m)')
    plt.xlim([0,L/m])
    plt.ylim(0,12)
    plt.text(35/m,10/(kN*m),r'$AF_{Moment}=\dfrac{%.03f \:\mathrm{kN-m}}{%.03f \:\mathrm{kN-m}} = %.03f$' % (MmaxP/(kN*m),Mmax1/(kN*m),MmaxP/Mmax1),fontsize=8)
    plt.savefig(os.path.join(save_folder, 'Figure_XX_OneWayExample_Moment'), dpi=300)

    # Shear Diagram
    plt.figure(figsize=(3.5,2.5))
    plt.axes(position=[0.15,0.17,0.82,0.80])
    plt.plot(resultsP.position_along_length/m,resultsP.shear_along_length/kN,label='With Ponding')
    plt.plot(results1.position_along_length/m,results1.shear_along_length/kN,label='No Ponding',color='tab:grey')
    plt.plot([0,0], [VmaxP/kN,Vmax1/kN],linestyle='None',marker='.',color='k')
    plt.plot([0,150/m],[VmaxP/kN,VmaxP/kN],'--',color='k')
    plt.plot([0,150/m],[Vmax1/kN,Vmax1/kN],'--',color='k')
    plt.text(95/m,(VmaxP-0.13)/kN,f'{VmaxP/kN:.03f}', fontsize=8)
    plt.text(95/m,(Vmax1-0.13)/kN,f'{Vmax1/kN:.03f}', fontsize=8)
    plt.legend()
    plt.xlabel('Position along length of beam (m)')
    plt.ylabel('Shear (kN)')
    plt.xlim([0,L/m])
    plt.ylim([-2,6])
    #plt.yticks(np.arange(-0.4,1.6,0.2))
    plt.text(200/m,0.5/kN,r'$AF_{Shear}=\dfrac{%.03f \:\mathrm{kN}}{%.03f \:\mathrm{kN}} = %.03f$' % (VmaxP/kN,Vmax1/kN,VmaxP/Vmax1),fontsize=8)
    plt.savefig(os.path.join(save_folder, 'Figure_XX_OneWayExample_Shear'), dpi=300)


plt.show()
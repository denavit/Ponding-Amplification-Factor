import numpy as np
import matplotlib.pyplot as plt
from math import pi, cos, cosh
from PyPonding.structures import ElasticBeam2d
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


# Input parameters
L       = 40*ft         # Beam span
S       = 5*ft          # Tributary width
E       = 29000.0*ksi   # Modulus of elasticity
I       = 75*inch**4    # Moment of inertia
gamma   = 62.4*pcf      # Unit weight of water
qD      = 0*psf         # Dead load
slope   = 0.25*in_per_ft   # Beam slope
C       = (gamma*S*L**4)/(pi**4*E*I)
Bpo     = 1/(1-C)


# Define loading cases
zw_over_zh_max = 2.0
n = 81 # number of analyses

zh = slope*L
zw = np.linspace(0,zw_over_zh_max*zh,n)

Bp_V  = np.zeros(n)
Bp_M  = np.zeros(n)
Bp_TL = np.zeros(n)


# Run Analyses
for i in range(n):
    # Build beam object
    beam = ElasticBeam2d(L,S,E,I,gamma,zj=zh,qD=qD)
    beam.num_elements = 160
    beam.maximum_number_of_iterations = 100
    beam.load_tolerance = 1e-6

    # Run Ponding Analysis
    resultsP = beam.run_analysis_OPS('IterativeLevel',target_zw=zw[i])
    VmaxP = np.amax(np.absolute(resultsP.shear_along_length))
    MmaxP = np.amax(resultsP.bending_moment_along_length)
    x_at_MmaxP = resultsP.position_along_length[np.argmax(resultsP.bending_moment_along_length)]
    TotalLoadP = resultsP.shear_along_length[0]-resultsP.shear_along_length[-1]    
    
    # Run First-Order Analysis
    beam.include_ponding_effect = False
    results1 = beam.run_analysis_OPS('IterativeLevel',target_zw=zw[i])
    Vmax1 = np.amax(np.absolute(results1.shear_along_length))
    Mmax1 = np.amax(results1.bending_moment_along_length)
    x_at_Mmax1 = results1.position_along_length[np.argmax(results1.bending_moment_along_length)]
    TotalLoad1 = results1.shear_along_length[0]-results1.shear_along_length[-1]
    
    Bp_V[i] = VmaxP/Vmax1
    Bp_M[i] = MmaxP/Mmax1
    Bp_TL[i] = TotalLoadP/TotalLoad1
   

# Plot Results
save_folder = 'figures'
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

plt.rc('axes',labelsize=8)
plt.rc('axes',titlesize=8)
plt.rc('legend',fontsize=8)
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)

plt.figure(figsize=(3.5,2.5))
plt.axes(position=[0.14,0.17,0.82,0.80])
#line0, = plt.plot([0.0,0.2,0.8,zw_over_zh_max],[1,1,Bpo,Bpo],'--', label='Idealized')
line0, = plt.plot([0.0,zw_over_zh_max],[Bpo,Bpo],'--', color='k', label='Basic')
line1, = plt.plot(zw/zh,Bp_M ,'-' ,label='$AF_{Moment}$')
line2, = plt.plot(zw/zh,Bp_TL,'-.',label='$AF_{Total\ Load}$')
line3, = plt.plot(zw/zh,Bp_V ,':' ,label='$AF_{Shear}$')
plt.xlabel('Normalized Water Level, $z_w/z_h$')
plt.ylabel('Amplification Factor')
plt.xlim(0.0,zw_over_zh_max)
plt.legend(handles=[line1,line2,line3,line0])
plt.savefig(os.path.join(save_folder, 'Figure_XX_OneWayExample_AF'), dpi=300)

plt.show()


# Print Data
print(f'C   = {C:.03f}')
print(f'Bpo = {Bpo:.03f}')
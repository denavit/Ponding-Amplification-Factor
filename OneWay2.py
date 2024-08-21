import os
import numpy as np
import matplotlib.pyplot as plt
from math import pi, cos, cosh
from PyPonding.structures import ElasticBeam2d

#This script performs an iterative ponding analysis for a beam with differing 
#support stiffness conditions at differing water levels.


# Plot Settings
plt.rc('axes',labelsize=8)
plt.rc('axes',titlesize=8)
plt.rc('legend',fontsize=8)
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)


# Define units
inch    = 1.0
kip     = 1.0
lb      = kip/1000.0
ft      = 12.0*inch
ksi     = kip/inch**2
psf     = lb/ft**2
pcf     = psf/ft
in_per_ft = inch/ft
mm      = inch/25.4
kN      = 224.80894*lb
m       = 1000*mm
kNm     = kN*m


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

# Define zw and zh for normalized water level
zh = slope*L
zw = np.linspace(0,zw_over_zh_max*zh,n)

# Define spring stiffnesses (kzi,kzj)
kz_list = [(None,None), (2,2), (None,2)]

# Initiate shear,moment, total load, x values, and amplification factor dicts
Bp_V_list  = dict()
Bp_M_list  = dict()
Bp_TL_list = dict()

x_values   = dict()
x_values1  = dict()

R1_list_AF = dict()
R2_list_AF = dict()

for j, kz in enumerate(kz_list):
    
    Bp_V  = np.zeros(n)
    Bp_M  = np.zeros(n)
    Bp_TL = np.zeros(n)
    
    R1_list      = np.zeros(n)
    R2_list      = np.zeros(n)
    x_resultant  = np.zeros(n)
    
    R1_list1      = np.zeros(n)
    R2_list1      = np.zeros(n)
    x_resultant1  = np.zeros(n)
    
    # Run Analyses
    for i in range(n):
        # Build beam object
        beam = ElasticBeam2d(L,S,E,I,gamma,zj=zh,qD=qD,kzi=kz[0],kzj=kz[1])
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
        
        Bp_V[i]  = VmaxP/Vmax1
        Bp_M[i]  = MmaxP/Mmax1
        Bp_TL[i] = TotalLoadP/TotalLoad1
        
        R1_list[i] = abs(resultsP.shear_along_length[0])
        R2_list[i] = abs(resultsP.shear_along_length[-1])
        x_resultant[i] = R2_list[i]*L/(R1_list[i]+R2_list[i])
        
        R1_list1[i] = abs(results1.shear_along_length[0])
        R2_list1[i] = abs(results1.shear_along_length[-1])
        x_resultant1[i] = R2_list1[i]*L/(R1_list1[i]+R2_list1[i])
        
    # Save values
    Bp_V_list[j]  = Bp_V
    Bp_M_list[j]  = Bp_M
    Bp_TL_list[j] = Bp_TL
    
    x_values[j]   = x_resultant
    x_values1[j]  = x_resultant1
    
    R2_list_AF[j] = R2_list/R2_list1
    R1_list_AF[j] = R1_list/R1_list1


# Plot Results
save_folder = os.path.join('figures', 'Paper_Figures')
if not os.path.exists(save_folder):
    os.makedirs(save_folder)


# Plot of amplfication factor for shear, moment, and total load vs zw/zh
# Also shows the basic amplficiation factor (if supports are rigid)
# One plot for each kz pair analyzed
for j, kz in enumerate(kz_list):
    plt.figure(figsize=(3.5,2.5))
    plt.axes(position=[0.16,0.17,0.80,0.75])
    line0, = plt.plot([0.0,zw_over_zh_max],[Bpo,Bpo],'--', color='k', label='Basic')
    line1, = plt.plot(zw/zh,Bp_V_list[j] ,':' ,label='$AF_{shear}$')
    line2, = plt.plot(zw/zh,Bp_M_list[j] ,'-' ,label='$AF_{moment}$')
    line3, = plt.plot(zw/zh,Bp_TL_list[j],'-.',label='$AF_{total\ load}$')
    plt.title(f'kzi = {kz[0]}, kzj = {kz[1]}')
    plt.xlabel('Normalized Water Level, $z_w/z_h$')
    plt.ylabel(f'Amplification Factor')
    plt.xlim(0.0,zw_over_zh_max)
    plt.legend(handles=[line1,line2,line3,line0])



# Normalized location of centroid wrt the left end vs zw/zh
plt.figure(figsize=(3.5,2.5))
plt.axes(position=[0.17,0.17,0.79,0.80])
line1, = plt.plot(zw/zh,x_values[0]/L,'-' ,label='Rigid Supports')
line2, = plt.plot(zw/zh,x_values[2]/L,'-.' ,label='Flexible Support on Right')
line3, = plt.plot(zw/zh,x_values1[0]/L,':' ,label='Without Ponding')
plt.xlabel('Normalized Water Level, $z_w/z_h$')
plt.ylabel('Location of Force Resultant\nfrom Left End Divided by Beam Length')
plt.xlim(0.0,zw_over_zh_max)
plt.ylim(0.0,0.5)
plt.legend()
plt.savefig(os.path.join(save_folder, 'Figure 09'), dpi=300)


# Right reaction AF vs zw/zh
plt.figure(figsize=(3.5,2.5))
plt.axes(position=[0.14,0.17,0.82,0.80])
line1, = plt.plot(zw/zh,R2_list_AF[0],'-' ,label='Rigid Supports')
line2, = plt.plot(zw/zh,R2_list_AF[1],'-.',label='Flexible Supports')
line3, = plt.plot(zw/zh,R2_list_AF[2],':' ,label='Flexible Support on Right')
plt.xlabel('Normalized Water Level, $z_w/z_h$')
plt.ylabel('Amplification of Right Support Reaction')
plt.xlim(0.0,zw_over_zh_max)
plt.legend()


#Left reaction AF vs zw/zh
plt.figure(figsize=(3.5,2.5))
plt.axes(position=[0.14,0.17,0.82,0.80])
line1, = plt.plot(zw/zh,R1_list_AF[0],'-' ,label='Rigid Supports')
line2, = plt.plot(zw/zh,R1_list_AF[1],'-.',label='Flexible Supports')
line3, = plt.plot(zw/zh,R1_list_AF[2],':' ,label='Flexible Support on Right')
plt.xlabel('Normalized Water Level, $z_w/z_h$')
plt.ylabel('Amplification of Left Support Reaction')
plt.xlim(0.0,zw_over_zh_max)
plt.legend()


# Right reaction with ponding and no ponding for rigid L and flex R case,
# also showing amplification factor.
plt.figure(figsize=(3.5,2.5))
ax1 = plt.axes(position=[0.16,0.17,0.70,0.80])
ax1.plot(zw/zh,R2_list/kN,label='With Ponding')
ax1.plot(zw/zh,R2_list1/kN,label='Without Ponding',color='tab:grey')
ax1.plot([],[],label='AF',linestyle=':',color='k')
ax1.legend()
ax1.set_xlabel('Normalized Water Level, $z_w/z_h$')
ax1.set_ylabel('Right Support Reaction (kN)')
ax1.set_xlim(0.0,zw_over_zh_max)
ax1.set_ylim(0,100)
ax2 = ax1.twinx()
ax2.set_ylabel('Amplification of Right Support Reaction')
ax2.plot(zw/zh,R2_list_AF[2], color='k',linestyle=':')
ax2.set_ylim(1,4)
plt.savefig(os.path.join(save_folder, 'Figure 10'), dpi=300)

plt.show()

# Print Data
print(f'C   = {C:.03f}')
print(f'Bpo = {Bpo:.03f}')
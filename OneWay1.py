import os
import numpy as np
import matplotlib.pyplot as plt
from math import pi, cos, cosh
from PyPonding.structures import ElasticBeam2d


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


def run_single_one_way_analysis(zw_over_zh, make_Figure_5=False):

    zh = slope*L
    zw = zw_over_zh*zh
    
    # Build beam object
    beam = ElasticBeam2d(L,S,E,I,gamma,zj=slope*L,qD=qD)
    beam.num_elements = 60
    
    # Run Ponding Analysis
    resultsP = beam.run_analysis_OPS('IterativeLevel',target_zw=zw)
    VmaxP = np.amax(np.absolute(resultsP.shear_along_length))
    MmaxP = np.amax(resultsP.bending_moment_along_length)
    x_at_MmaxP = resultsP.position_along_length[np.argmax(resultsP.bending_moment_along_length)]
    TotalLoadP = resultsP.shear_along_length[0]-resultsP.shear_along_length[-1]
    
    #beam.plot_deformed(zw=zw)
    
    # Run First-Order Analysis
    beam.include_ponding_effect = False
    results1 = beam.run_analysis_OPS('IterativeLevel',target_zw=zw)
    Vmax1 = np.amax(np.absolute(results1.shear_along_length))
    Mmax1 = np.amax(results1.bending_moment_along_length)
    x_at_Mmax1 = results1.position_along_length[np.argmax(results1.bending_moment_along_length)]
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
   
    AF = max(MmaxP/Mmax1,VmaxP/Vmax1)
    
    if make_Figure_5:
        
        plt.figure(figsize=(3.5,2.5))
        plt.axes(position=(0.14,0.17,0.84,0.80))
        plt.plot(resultsP.position_along_length/m,resultsP.shear_along_length/kN,label='With Ponding')
        plt.plot(results1.position_along_length/m,results1.shear_along_length/kN,label='Without Ponding',color='tab:grey')
        plt.plot([0,0], [VmaxP/kN,Vmax1/kN],linestyle='None',marker='.',color='k')
        plt.plot([0,0.20*L/m],[VmaxP/kN,VmaxP/kN],'--',color='k')
        plt.plot([0,0.20*L/m],[Vmax1/kN,Vmax1/kN],'--',color='k')
        plt.text(0.20*L/m,VmaxP/kN,f' {VmaxP/kN:.3f} kN', va='center', fontsize=8)
        plt.text(0.20*L/m,Vmax1/kN,f' {Vmax1/kN:.3f} kN', va='center', fontsize=8)
        plt.legend(loc='upper right')
        plt.xlabel('Position along length of beam (m)')
        plt.ylabel('Shear (kN)')
        plt.xlim([0,L/m])
        plt.ylim([-2,6.5])
        plt.text(5,0.5,r'$AF_{shear}=\dfrac{'+f'{VmaxP/kN:.3f}'+r'\text{ kN}}{'+f'{Vmax1/kN:.3f}'+r'\text{ kN}}='+f'{VmaxP/Vmax1:.3f}'+r'$',fontsize=8)
        plt.savefig(os.path.join(save_folder, 'Figure 05a - Shear'), dpi=300)        
        
        plt.figure(figsize=(3.5,2.5))
        plt.axes(position=(0.14,0.17,0.84,0.80))
        plt.plot(resultsP.position_along_length/m,resultsP.bending_moment_along_length/kNm,label='With Ponding')
        plt.plot(results1.position_along_length/m,results1.bending_moment_along_length/kNm,label='Without Ponding',color='tab:grey')
        plt.plot([x_at_MmaxP/m,x_at_Mmax1/m],[MmaxP/kNm,Mmax1/kNm],linestyle='None',marker='.',color='k')
        plt.plot([0,x_at_MmaxP[0]/m],[MmaxP/kNm,MmaxP/kNm],'--',color='k')
        plt.plot([0,x_at_Mmax1[0]/m],[Mmax1/kNm,Mmax1/kNm],'--',color='k')
        plt.text(0.1,MmaxP/kNm+0.25,f' {MmaxP/kNm:.3f} kN-m', fontsize=8)
        plt.text(2.2,Mmax1/kNm+0.25,f' {Mmax1/kNm:.3f} kN-m', fontsize=8)
        plt.legend(loc='upper right')
        plt.xlabel('Position along length of beam (m)')
        plt.ylabel('Bending Moment (kN-m)')
        plt.xlim([0,L/m])
        plt.ylim([0,12])
        plt.text(0.9,1.8,r'$AF_{moment}=\dfrac{'+f'{MmaxP/kNm:.3f}'+r'\text{ kN-m}}{'+f'{Mmax1/kNm:.3f}'+r'\text{ kN-m}}='+f'{MmaxP/Mmax1:.3f}'+r'$',fontsize=8)
        plt.savefig(os.path.join(save_folder, 'Figure 05b - Moment'), dpi=300)
        
        plt.show()
        
    return AF
        
def Figure_6():
    
    zw_over_zh_max = 2.0
    zw_over_zh_list = np.linspace(0.001,zw_over_zh_max,100)
    
    Cs = (gamma*S*L**4)/(pi**4*E*I)
    basic_AF = 1/(1-Cs)
    
    AF_list = []
    for zw_over_zh in zw_over_zh_list:
        AF = run_single_one_way_analysis(zw_over_zh)
        AF_list.append(AF)

    plt.figure(figsize=(3.5,2.5))
    plt.axes(position=(0.14,0.17,0.82,0.80))
    plt.plot(zw_over_zh_list,AF_list,color='tab:grey', label='Iterative Analysis')
    plt.plot([0,zw_over_zh_max],[basic_AF,basic_AF],color='k', linestyle='--',label='Basic Amplification')
    plt.plot([0,0.15,0.75,zw_over_zh_max],[1,1,basic_AF,basic_AF],color='r', linestyle=':',label='Considering Water Height')   
    plt.legend()
    plt.xlabel('Normalized Water Height, $z_w/z_h$')
    plt.ylabel('Amplification Factor')
    plt.xlim([0,zw_over_zh_max])
    plt.ylim([0.9,2.3])
    plt.savefig(os.path.join(save_folder, 'Figure 06'), dpi=300, bbox_inches='tight')

    plt.show()


if __name__ == "__main__":
    #This script will perform an iterative ponding analysis for a simply supported
    #sloped beam for many water levels. Sample single analysis plots for shear and 
    #moment values are generated to illustrate the drivation of the amplification 
    #factor. The change in amplification for many water levels is also plotted.

    save_folder = os.path.join('figures', 'Paper_Figures')
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    # Input parameters
    L       = 40*ft          # Beam span
    S       = 5*ft           # Tributary width
    E       = 29000.0*ksi    # Modulus of elasticity
    I       = 75*inch**4     # Moment of inertia
    gamma   = 62.4*pcf       # Unit weight of water
    slope   = 0.25*in_per_ft # Beam slope
    qD      = 0*psf          # Dead load

    # Run one analysis with zw/zh = 0.5 to make Figure 5
    run_single_one_way_analysis(zw_over_zh=0.5,make_Figure_5=True)

    # Run multiple analyses with varying zw/zh to make Figure 6
    Figure_6()

import matplotlib.pyplot as plt
import os
from TwoWay_RunAnalysis import run_single_analysis

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


# Run Analyses
zw_over_zh_list = [0.001,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
                         0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00,
                         1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,1.50,
                         1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,1.95,2.00]
#zw_over_zh_list = [0.001,0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00]

MmaxP_list = []
Mmax0_list = []
MmaxP_over_Mmax0_list = []
for zw_over_zh in zw_over_zh_list:
    (amplification_factors,extra_output) = run_single_analysis(case='A',Cp=0.3,Cs=0.3,roof_slope=0.5*in_per_ft,zw_over_zh=zw_over_zh,num_spaces=16,qD=0,return_extra_output=True)
    MmaxP_list.append(extra_output['MmaxP']/kN/1000)
    Mmax0_list.append(extra_output['Mmax0']/kN/1000)
    MmaxP_over_Mmax0_list.append(extra_output['MmaxP']/extra_output['Mmax0'])


# Plot Results
save_folder = os.path.join('figures', 'Paper_Figures')
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

plt.figure(figsize=(3.5,2.5))
ax1 = plt.axes(position=[0.15,0.17,0.71,0.78])
ax1.plot(zw_over_zh_list,MmaxP_list,label='With Ponding')
ax1.plot(zw_over_zh_list,Mmax0_list,label='Without Ponding',color='tab:grey')
ax1.plot([],[],label='Amplification',linestyle=':',color='k')
ax1.legend(frameon=False)
ax1.set_xlabel('Normalized Water Level, $z_w/z_h$')
ax1.set_ylabel(r'Bending Moment (kN-m$\times 10^3$)')
ax1.set_xlim(0.0,max(zw_over_zh_list))
ax1.set_ylim(0,10)
ax2 = ax1.twinx()
ax2.set_ylabel('Amplification of Bending Moment')
ax2.plot(zw_over_zh_list,MmaxP_over_Mmax0_list, color='k',linestyle=':')
ax2.set_ylim(1,2.5)
plt.savefig(os.path.join(save_folder, 'Figure 12'), dpi=300)

plt.show()

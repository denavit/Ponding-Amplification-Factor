from PyPonding.structures import IdealizedBay
from math import pi
import json

# Define units
inch      = 1.0
kip       = 1.0
lb        = kip/1000.0
ft        = 12.0*inch
in_per_ft = inch/ft
ksi       = kip/inch**2
plf       = lb/ft
psf       = lb/ft**2
pcf       = psf/ft
kipft     = kip*ft

# Define members to investigate for each case
titles = dict()
titles['Flat']  = ["Primary Members", "Secondary Members", "Total Load"]
titles['A']     = ["Top Primary Member", "Secondary Members"]
titles['B']     = ["Top Primary Member", "Bottom Primary Member", "Secondary Members"]
titles['C']     = ["Primary Members", "Secondary Member 2"]
titles['D']     = ["Primary Members", "Secondary Member 1"]
titles['E']     = ["Top Primary Member", "Bottom Primary Member", "Secondary Member 1"]
titles['F']     = ["Bottom Primary Member", "Secondary Member 2"]


def run_single_analysis(case,Cp,Cs,roof_slope,zw_over_zh,qD):

    # Bay Geometry
    Lp = 40*ft
    Ls = 40*ft
    num_spaces = 16
    S  = Lp/num_spaces

    # Loads
    γw = 62.4*pcf

    # Material properties
    E = 29000.0*ksi

    # Member stiffness
    Ip = (γw*Ls*Lp**4)/(pi**4*E*Cp)
    Is = (γw*S*Ls**4)/(pi**4*E*Cs)
    
    if case == 'Flat':
        
        zw = zw_over_zh # Special case without slope 
        
        z_TL = 0.0*inch
        z_TR = 0.0*inch
        z_BL = 0.0*inch
        z_BR = 0.0*inch
        
        edge_condition_L = 'mirrored'
        edge_condition_R = 'mirrored'
        edge_condition_T = 'mirrored'
        edge_condition_B = 'mirrored'

        
    elif case == 'A':
        
        zh = roof_slope*Lp
        zw = zw_over_zh*zh
        
        z_TL = zh
        z_TR = zh
        z_BL = 0.0*inch
        z_BR = 0.0*inch
        
        edge_condition_L = 'mirrored'
        edge_condition_R = 'mirrored'
        edge_condition_T = 'mirrored'
        edge_condition_B = 'rigid'
        
    elif case == 'B':
        
        zh = roof_slope*Lp
        zw = zw_over_zh*zh
        
        z_TL = 0.0*inch
        z_TR = 0.0*inch
        z_BL = zh
        z_BR = zh
        
        edge_condition_L = 'mirrored'
        edge_condition_R = 'mirrored'
        edge_condition_T = 'mirrored'
        edge_condition_B = 'mirrored'
        
    elif case == 'C':
        
        zh = roof_slope*Ls
        zw = zw_over_zh*zh
        
        z_TL = 0.0*inch
        z_TR = zh
        z_BL = 0.0*inch
        z_BR = zh
        
        edge_condition_L = 'rigid'
        edge_condition_R = 'mirrored'
        edge_condition_T = 'mirrored'
        edge_condition_B = 'mirrored'
        
    elif case == 'D':
        
        zh = roof_slope*Ls
        zw = zw_over_zh*zh
        
        z_TL = 0.0*inch
        z_TR = zh
        z_BL = 0.0*inch
        z_BR = zh
        
        edge_condition_L = 'mirrored'
        edge_condition_R = 'mirrored'
        edge_condition_T = 'mirrored'
        edge_condition_B = 'mirrored'
        
    elif case == 'E':
        
        zh = roof_slope*Lp
        zw = zw_over_zh*zh
        
        z_TL = 0.0*inch
        z_TR = zh
        z_BL = zh
        z_BR = zh
        
        edge_condition_L = 'mirrored'
        edge_condition_R = 'mirrored'
        edge_condition_T = 'mirrored'
        edge_condition_B = 'mirrored'
        
    elif case == 'F':
        
        zh = roof_slope*Lp
        zw = zw_over_zh*zh
        
        z_TL = 0.0*inch
        z_TR = 0.0*inch
        z_BL = 0.0*inch
        z_BR = zh
        
        edge_condition_L = 'rigid'
        edge_condition_R = 'mirrored'
        edge_condition_T = 'rigid'
        edge_condition_B = 'mirrored'
        
    else:
        raise ValueError(f'Unknown case {case}')

    # Define IdealizedBay object
    bay_input = {
        'primary_member_span': Lp,
        'secondary_member_span': Ls,
        'number_of_joist_spaces': num_spaces,
        'dead_load_uniform': qD,
        'dead_load_on_top_primary_member': 0.0*plf,
        'dead_load_on_bottom_primary_member': 0.0*plf,
        'water_density': γw,  
        'snow_density': 0.0*pcf,
        'snow_height': 0.0*inch,
        'alpha': 1.0,
        'load_factor_dead':    1.0,
        'load_factor_ponding': 1.0,
        'load_factor_snow':    1.0,
        'z_TL': z_TL,
        'z_TR': z_TR,
        'z_BL': z_BL,
        'z_BR': z_BR,
        'secondary_member_camber': 0.000*inch,
        'primary_member_camber_T': 0.000*inch,
        'primary_member_camber_B': 0.000*inch,
        'edge_condition_L': edge_condition_L,
        'edge_condition_R': edge_condition_R,
        'edge_condition_T': edge_condition_T,
        'edge_condition_B': edge_condition_B,
        'E':  E,
        'As': 100*inch**2,
        'Ap': 100*inch**2,
        'Is': Is,
        'Ip': Ip,
        'analsis_engine': 'FE',
    }
    bay = IdealizedBay(**bay_input)

    bay.include_ponding_effect = False
    results0 = bay.Run_Analysis(zw)
    bay.include_ponding_effect = True
    resultsP = bay.Run_Analysis(zw)

    amplification_factors = dict()

    AF = 1/(1-1.15*Cp-Cs)
    print('\nBasic Amplification 1/(1-1.15Cp-Cs) = %.3f' % (AF))
    amplification_factors['basic'] = AF
    
    print('==== Total Load ====')
    print(f'No Ponding Effect:    {results0.total_factored_load:.3f}')
    print(f'With Ponding Effect:  {resultsP.total_factored_load:.3f}')
    print(f'Amplification Factor: {resultsP.total_factored_load/results0.total_factored_load:.3f}')
    amplification_factors['total_load'] = resultsP.total_factored_load/results0.total_factored_load

    print('==== Member Output ====')
    print('  Member         Max Shear      Max Moment         Max Shear      Max Moment        Shear AF         Moment AF')
    print('                  (kips)         (kip-ft)           (kips)         (kip-ft)           (---)            (---)  ')
    
    # Secondary Members
    secondary_amplification_factors = []
    for i in range(bay.number_of_joist_spaces+1):
        if i == 0 and bay.edge_condition_L == 'rigid':
            secondary_amplification_factors.append(None)
            continue
        if i == bay.number_of_joist_spaces and bay.edge_condition_R == 'rigid':
            secondary_amplification_factors.append(None)
            continue  
        Vmax0 = max(results0.secondary_members_bot_reaction[i],results0.secondary_members_top_reaction[i])
        Mmax0 = max(results0.secondary_members_moment[i,:])/12
        VmaxP = max(resultsP.secondary_members_bot_reaction[i],resultsP.secondary_members_top_reaction[i])
        MmaxP = max(resultsP.secondary_members_moment[i,:])/12
        if Vmax0 == 0.0:
            if VmaxP == 0.0:
                V_AF_str = '  Undef.'
            else:
                V_AF_str = '   Inf  '
        else:
            V_AF_str = f'{VmaxP/Vmax0:>8.3f}'
        if Mmax0 == 0.0:
            if MmaxP == 0.0:
                M_AF_str = '  Undef.'
                secondary_amplification_factors.append('Undefined')
            else:
                M_AF_str = '   Inf  '
                secondary_amplification_factors.append('Infinite')
        else:
            M_AF_str = f'{MmaxP/Mmax0:>8.3f}'
            secondary_amplification_factors.append(max(VmaxP/Vmax0,MmaxP/Mmax0))
        
        print(f'Secondary {i+1:<2d}   {Vmax0:>8.2f}         {Mmax0:>8.2f}         {VmaxP:>8.2f}         {MmaxP:>8.2f}         {V_AF_str}         {M_AF_str}')
    
    amplification_factors['secondary'] = secondary_amplification_factors
        
        
    # Top Primary Member
    if bay.edge_condition_T == 'rigid':
        amplification_factors['top'] = None
    else:
        Vmax0 = max(abs(results0.top_primary_member_shear))
        Mmax0 = max(results0.top_primary_member_moment)/12
        VmaxP = max(abs(resultsP.top_primary_member_shear))
        MmaxP = max(resultsP.top_primary_member_moment)/12
        print(f'Primary Top    {Vmax0:>8.2f}         {Mmax0:>8.2f}         {VmaxP:>8.2f}         {MmaxP:>8.2f}         {VmaxP/Vmax0:>8.3f}         {MmaxP/Mmax0:>8.3f}')
        amplification_factors['top'] = max(VmaxP/Vmax0,MmaxP/Mmax0)
        
    # Bottom Primary Member
    if bay.edge_condition_B == 'rigid':
        amplification_factors['bottom'] = None
    else:
        Vmax0 = max(abs(results0.bot_primary_member_shear))
        Mmax0 = max(results0.bot_primary_member_moment)/12
        VmaxP = max(abs(resultsP.bot_primary_member_shear))
        MmaxP = max(resultsP.bot_primary_member_moment)/12
        print(f'Primary Bottom {Vmax0:>8.2f}         {Mmax0:>8.2f}         {VmaxP:>8.2f}         {MmaxP:>8.2f}         {VmaxP/Vmax0:>8.3f}         {MmaxP/Mmax0:>8.3f}')
        amplification_factors['bottom'] = max(VmaxP/Vmax0,MmaxP/Mmax0)
  
    return amplification_factors

def run_analysis_loop(case,Cp_list,Cs_list,roof_slope,zw_over_zh_list,qD):
    '''run multiple analyses and save results to json'''
    
    # Get percentage completion
    count = 0
    total = len(Cs_list) * len(Cp_list) * len(zw_over_zh_list)
    
    # Initalize Bp list
    Bp_results = dict()
    
    # Loop Cs values
    for Cs in Cs_list:
                 
        # Initalize Bp list
        Bp_results[Cs] = dict()
     
        # Run analyses for each Cp value
        for iCp, Cp in enumerate(Cp_list):
            
            # Inititalize results data structure
            Bp_results[Cs][Cp] = dict()
            for title in titles[case]:
                Bp_results[Cs][Cp][title] = []

            # Run analyises
            for zw_over_zh in zw_over_zh_list:
                
                # Get percentage completion
                count += 1
                completion_status = round(100 * (count/total), 1)
                
                print(f'\n==== Analysis {completion_status}% complete. ==== \n==== Analysis for Case = {case}, {Cs = }, {Cp = }, {zw_over_zh = } ====\n')
                
                amplification_factors = run_single_analysis(case,Cp,Cs,roof_slope,zw_over_zh,qD)       
                                                                      
                # Save needed results
                for title in titles[case]:
                    if title == "Primary Members":
                        Bp_max = max(amplification_factors['top'], amplification_factors['bottom'])
                    elif title == "Top Primary Member":
                        Bp_max = amplification_factors['top']
                    elif title == "Bottom Primary Member":
                        Bp_max = amplification_factors['bottom']
                    elif title == "Secondary Members":
                        Bp_max = max(amplification_factors['secondary'])
                    elif title == "Secondary Member 1":
                        Bp_max = amplification_factors['secondary'][0]
                    elif title == "Secondary Member 2":
                        Bp_max = amplification_factors['secondary'][1]
                    elif title == "Total Load":
                        Bp_max = amplification_factors['total_load']
                    else:
                        return print(f"Unknown title: '{title}'") 
                    
                    Bp_results[Cs][Cp][title].append(Bp_max) 
    
    with open(f'Bp_results_{case}.json', 'w') as j:
    # write the dictionary to the file in JSON format
        json.dump(Bp_results, j)
                
    return Bp_results


# Run Analyses
roof_slope = 0.5*in_per_ft
qD = 0.0*psf
Cs_list = [0.001, 0.1, 0.2, 0.3]
Cp_list = [0.001, 0.1, 0.2, 0.3]
zw_over_zh_list = [0.001,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
                    0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,1.00,1.05,1.10,
                    1.15,1.20,1.25,1.30,1.35,1.40,1.45,1.50]

cases = ['Flat','A','B','C','D','E','F']
for case in cases:
    Bp_results = run_analysis_loop(case,Cp_list,Cs_list,roof_slope,zw_over_zh_list,qD)
    print(Bp_results)

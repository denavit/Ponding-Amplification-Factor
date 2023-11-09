import json
from pprint import pprint

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


def equation_generator(case,Cs_list,Cp_list,Cs_coefs=None,Cp_coefs=None): 
    '''Read resutls generated by TwoWay_RunAnalysis.py and determine the optimal coefficients'''
    
    Bp_results = json.load(open(f'Bp_results_{case}.json'))
    
    if Cs_coefs is None:
        # Set to default value
        Cs_coefs = [0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
                    0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00,
                    1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,1.50,
                    1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,1.95,2.00]
    
    if Cp_coefs is None:
        # Set to default value
        Cp_coefs = [0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
                    0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00,
                    1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,1.50,
                    1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,1.95,2.00]
        
    # Initialize dictionary of coefficients
    ideal_coefs = dict()
    
    # Find equation for each member category
    for title in titles[case]:
        
        ideal_coefs[title] = dict()  
        
        error_sum = float('inf') # inf so that first error automatically replaces it    
        for Cs_coef in Cs_coefs:
            for Cp_coef in Cp_coefs:
                                   
                # Compute error for this combination of coefs
                test_error_sum = 0
                for Cs in Cs_list:                    
                    for Cp in Cp_list:
                        Bp_ref = max(Bp_results[str(Cs)][str(Cp)][title]) # Takes the max amplification
                        
                        if (1-Cp_coef*Cp-Cs_coef*Cs) <= 0:
                            error = float('inf')
                        else:
                            M = 1
                            if Bp_ref-(1/(1-Cp_coef*Cp-Cs_coef*Cs)) > 0:
                                M = 5 # Penalty for underestimating the amplification factor
                            error = M*(Bp_ref-(1/(1-Cp_coef*Cp-Cs_coef*Cs)))**2
                        
                        test_error_sum += error
                        
                # Set this combination of coefs if it gives the lowest error
                if test_error_sum < error_sum:
                    ideal_coefs[title]['Cp_coef'] = Cp_coef
                    ideal_coefs[title]['Cs_coef'] = Cs_coef
                    error_sum = test_error_sum
    
    return ideal_coefs


# Run equation generator
Cs_list = [0.001, 0.1, 0.2, 0.3]
Cp_list = [0.001, 0.1, 0.2, 0.3]

cases = ['Flat','A','B','C','D','E','F']
for case in cases:
    ideal_coefs = equation_generator(case,Cs_list,Cp_list,Cs_coefs=None,Cp_coefs=None)
    print(f'\nCase {case}')
    pprint(ideal_coefs)

# Run this after the results files are created by TwoWay_RunAnalysis.py
# to generate the ideal coefficients for the amplificaiton factor equation

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


def equation_generator(case, num_spaces, member_title, Cs_coefs=None,Cp_coefs=None,target_zw_over_zh=None): 

    # Read results from json file
    results = json.load(open(f'results_{case}.json'))
    
    Cp_list = []
    Cs_list = []
    amplification_factor_list = dict()
    zw_over_zh_list = dict()
    
    for result in results:
        if result["num_spaces"] == num_spaces and member_title in result["amplification_factor"]:
            Cs_list.append(result["Cs"])
            Cp_list.append(result["Cp"])
            amplification_factor_list[f'Cs={result["Cs"]},Cp={result["Cp"]}'] = result["amplification_factor"][member_title]
            zw_over_zh_list[f'Cs={result["Cs"]},Cp={result["Cp"]}'] = result["zw_over_zh_list"]
    
    if Cs_coefs is None:
        # Set to default value
        Cs_coefs = [0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
                    0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00,
                    1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,1.50,
                    1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,1.95,2.00]

        if member_title in ["Secondary Members","Secondary Member 1","Secondary Member 2"]:
            Cs_coefs = [1.00,1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,1.50,
                             1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,1.95,2.00]

    if Cp_coefs is None:
        # Set to default value
        Cp_coefs = [0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
                    0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00,
                    1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,1.50,
                    1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,1.95,2.00]
        
    # initialize dictionary of coefficients
    ideal_coefs = dict()
    
    # find equation coefs
    
    error_sum = float('inf') # inf so that first error automatically replaces it    
    for Cs_coef in Cs_coefs:
        for Cp_coef in Cp_coefs:
                               
            # Compute error for this combination of coefs
            test_error_sum = 0
            for Csi, Cs in enumerate(Cs_list):                    
                for Cpi, Cp in enumerate(Cp_list):
                
                    # Determine reference amplification value
                    if target_zw_over_zh is None:
                        # Use the maximum amplification over all zw/zh
                        Bp_ref = max(amplification_factor_list[f'Cs={Cs},Cp={Cp}'])
                    else:
                        # Use the amplification for the target zw/zh
                        ind = zw_over_zh_list[f'Cs={Cs},Cp={Cp}'].index(target_zw_over_zh)
                        Bp_ref = amplification_factor_list[f'Cs={Cs},Cp={Cp}'][ind]
                    
                    # Calculate the incremenet in total error from this data point
                    if (1-Cp_coef*Cp-Cs_coef*Cs) <= 0:
                        error = float('inf')
                    else:
                        M = 1
                        if Bp_ref-(1/(1-Cp_coef*Cp-Cs_coef*Cs)) > 0:
                            M = 5 # Penalty for underestimating the amplification factor
                        error = M*(Bp_ref-(1/(1-Cp_coef*Cp-Cs_coef*Cs)))**2
                    
                    # Update total error
                    test_error_sum += error
                    
            # Set this combination of coefs if it gives the lowest error
            if test_error_sum < error_sum:
                ideal_coefs['Cp_coef'] = Cp_coef
                ideal_coefs['Cs_coef'] = Cs_coef
                error_sum = test_error_sum

    return ideal_coefs


def gen_all_equations(case_list=None): 
    #generates all equations for all members in all given cases for all member spacings
    
    if case_list is None:
        # Default to run all cases
        case_list = ["Flat","A","B","C","D","E","F"]
    
    # Convert case_list to list if necessary
    if type(case_list) is not list:
        case_list = [case_list]
    
    # Define members to investigate for each case
    titles = dict()
    titles['Flat']  = ["Primary Members", "Secondary Members", "Total Load"]
    titles['A']     = ["Top Primary Member", "Secondary Members"]
    titles['B']     = ["Top Primary Member", "Bottom Primary Member", "Secondary Members"]
    titles['C']     = ["Primary Members", "Secondary Member 2"]
    titles['D']     = ["Primary Members", "Secondary Member 1"]
    titles['E']     = ["Top Primary Member", "Bottom Primary Member", "Secondary Member 1"]
    titles['F']     = ["Top Primary Member", "Secondary Member 2"]
    
    # Find number of spaces
    num_spaces_list = [2,16]
    
    # Run equation generator and print results to screen
    for case in case_list:
        print(f'\n###########################\n### Case: {case:4}          ###\n###########################')
        for member_title in titles[case]:
            print(f'\n Member: {member_title}')
            print(' --------------------------')
            for num_spaces in num_spaces_list:
                ideal_coefs = equation_generator(case, num_spaces, member_title, Cs_coefs=None,Cp_coefs=None)
                print(f'  {num_spaces = :2}, xp = {ideal_coefs["Cp_coef"]:.2f}, xs = {ideal_coefs["Cs_coef"]:.2f}')


if __name__ == "__main__":

    '''                
    # Generate single set of coefficients           
    case = 'Flat'
    num_spaces = 16
    member_title = 'Total Load'
    ideal_coefs = equation_generator(case, num_spaces, member_title)
    print(ideal_coefs)
    '''
    
    '''
    ideal_coefs = equation_generator('B', 16, 'Top Primary Member')
    print(ideal_coefs)
    ideal_coefs = equation_generator('B', 16, 'Top Primary Member', target_zw_over_zh=1.5)
    print(ideal_coefs)
    '''
    
    gen_all_equations()





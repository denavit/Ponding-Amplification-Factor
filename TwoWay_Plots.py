# After running TwoWayAnalysis.py, this code will use the exported json file to
# create plots of the amplification factor vs the normalize water level for many 
# flexibility coefficients and all member types of a given bay configuration.

import matplotlib.pyplot as plt
import os
import json


# Plot formating options
plt.rc('font',family='serif')
plt.rc('mathtext',fontset='dejavuserif')
plt.rc('axes',labelsize=8)
plt.rc('axes',titlesize=8)
plt.rc('legend',fontsize=6)
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)
  

# Amplification factor rules dictionary 
Bp_rules = dict()
Bp_rules[('Flat','Primary Members')]    = {'Bpmin':None,  'zwzha':None, 'zwzhb':0.00, 'xp':1.00, 'xs':0.85}
Bp_rules[('Flat','Secondary Members')]  = {'Bpmin':None,  'zwzha':None, 'zwzhb':0.00, 'xp':1.15, 'xs':1.00}
Bp_rules[('Flat','Total Load')]         = {'Bpmin':None,  'zwzha':None, 'zwzhb':0.00, 'xp':0.85, 'xs':0.85}
Bp_rules[('A','Top Primary Member')]    = {'Bpmin':1.00,  'zwzha':0.20, 'zwzhb':0.85, 'xp':0.60, 'xs':1.15}
Bp_rules[('A','Secondary Members')]     = {'Bpmin':1.00,  'zwzha':0.20, 'zwzhb':0.85, 'xp':0.55, 'xs':1.00}
Bp_rules[('B','Top Primary Member')]    = {'Bpmin':1.00,  'zwzha':0.00, 'zwzhb':0.85, 'xp':1.20, 'xs':1.10}
Bp_rules[('B','Bottom Primary Member')] = {'Bpmin':1.00,  'zwzha':0.00, 'zwzhb':0.85, 'xp':1.00, 'xs':0.85}
Bp_rules[('B','Secondary Members')]     = {'Bpmin':1.00,  'zwzha':0.00, 'zwzhb':0.85, 'xp':1.15, 'xs':1.00}
Bp_rules[('C','Primary Members')]       = {'Bpmin':'Bp0', 'zwzha':0.00, 'zwzhb':0.85, 'xp':1.00, 'xs':0.85}
Bp_rules[('C','Secondary Member 2')]    = {'Bpmin':None,  'zwzha':None, 'zwzhb':0.00, 'xp':0.55, 'xs':1.00}
Bp_rules[('D','Primary Members')]       = {'Bpmin':'Bp0', 'zwzha':0.00, 'zwzhb':0.85, 'xp':1.00, 'xs':0.85}
Bp_rules[('D','Secondary Member 1')]    = {'Bpmin':None,  'zwzha':None, 'zwzhb':0.00, 'xp':0.30, 'xs':1.00}    
Bp_rules[('E','Top Primary Member')]    = {'Bpmin':1.00,  'zwzha':0.20, 'zwzhb':0.85, 'xp':1.15, 'xs':1.10}
Bp_rules[('E','Bottom Primary Member')] = {'Bpmin':1.00,  'zwzha':0.20, 'zwzhb':0.85, 'xp':1.00, 'xs':0.85}
Bp_rules[('E','Secondary Member 1')]    = {'Bpmin':None,  'zwzha':None, 'zwzhb':0.00, 'xp':0.30, 'xs':1.00}
Bp_rules[('F','Top Primary Member')]    = {'Bpmin':'Bp0', 'zwzha':0.00, 'zwzhb':0.60, 'xp':0.65, 'xs':0.95}
Bp_rules[('F','Secondary Member 2')]    = {'Bpmin':None,  'zwzha':None, 'zwzhb':0.00, 'xp':0.30, 'xs':1.00}


def two_way_plot(case, Cs, num_spaces, member_title, save_name=None, save_folder_name=None):
    # Make a plot of amplifcation factor vs zw/zh from results in json file

    zwzh_max = 1.5
    color_list = ['tab:blue','tab:orange','tab:green','tab:red']

    nice_member_title = {'Primary Members':'Primary',
                         'Secondary Members':'Secondary',
                         'Total Load':'Total Load',
                         'Top Primary Member':'Higher Primary',
                         'Bottom Primary Member':'Lower Primary',
                         'Secondary Member 1':'Lowest Secondary',
                         'Secondary Member 2':'Lowest Secondary'}


    # Create folder to save figures to if it doesn't exist
    if save_folder_name is None:
        save_folder_name = f'Config {case} - {nice_member_title[member_title]}'
    
    try:
        os.mkdir(os.path.join('figures',save_folder_name))
    except FileExistsError:
        pass    


    # Load json results
    results = json.load(open(f'results_{case}.json'))
    
    # Create figure
    fig = plt.figure(figsize=(3.25,2.75))
    ax = fig.add_axes([0.18,0.22,0.76,0.74])
   
    i = 0
    for result in results:
        if result["Cs"] != Cs or result["num_spaces"] != num_spaces:
            continue 
        
        Cp = result["Cp"]
        
        # Plot computed results
        x = result["zw_over_zh_list"]
        y = result["amplification_factor"][member_title]
        plt.plot(x,y, label=f'$C_p$ = {Cp}', color=color_list[i])  
        
        # Plot equation
        rules = Bp_rules[(case,member_title)]
        Bp_max = 1/(1-rules['xp']*Cp-rules['xs']*Cs)
        
        if rules['Bpmin'] is None:
            x = [0,zwzh_max]
            y = [Bp_max,Bp_max]
        else:
            if rules['Bpmin'] == 'Bp0':
                Bp_min = 1/(1-rules['xs']*Cs)
            else:
                Bp_min = rules['Bpmin']
            x = [0,rules['zwzha'],rules['zwzhb'],zwzh_max]
            y = [Bp_min,Bp_min,Bp_max,Bp_max]
        
        plt.plot(x,y,'--', color=color_list[i])
        
        i += 1
        
    # Plot formatting
    (_,ytop) = ax.get_ylim()
    if ytop <= 1.25:
        ax.set_yticks([1.00,1.05,1.10,1.15,1.20,1.25])
    #plt.title(f'Config. {case} --- $C_s$ = {Cs} --- {member_title} --- {num_spaces} spaces', fontsize = 8)
    plt.xlim((0,zwzh_max))
    plt.ylim(bottom=1,top=ytop)
    plt.xlabel('Normalized Water Level, $z_w/z_h$')
    plt.ylabel('Amplification')
    plt.legend(loc='upper center',frameon=True,ncol=4,bbox_to_anchor=(0.42,-0.18))

    # Save Figure
    if save_name is None:
        save_name = f'Cs_{Cs}_{num_spaces}spaces.png'
        
    plt.savefig(os.path.join('figures',save_folder_name, save_name), dpi = 300)

    # Close figure
    plt.close(fig)

def create_paper_figures():
    # Creates figures that appear in the paper named by the figure number.
    two_way_plot('Flat', 0.3, 16, 'Primary Members',    'Figure 08',  'Paper_Figures')
    two_way_plot(   'A', 0.3, 16, 'Secondary Members',  'Figure 12',  'Paper_Figures')
    two_way_plot(   'A', 0.3, 16, 'Top Primary Member', 'Figure 13',  'Paper_Figures')
    two_way_plot(   'B', 0.3, 16, 'Secondary Members',  'Figure 14',  'Paper_Figures')
    two_way_plot(   'C', 0.3, 16, 'Primary Members',    'Figure 15',  'Paper_Figures')
    two_way_plot(   'C', 0.3,  2, 'Secondary Member 2', 'Figure 16a', 'Paper_Figures')
    two_way_plot(   'C', 0.3, 16, 'Secondary Member 2', 'Figure 16b', 'Paper_Figures')
    
def create_all_figures(case_list=None):
    # Creates figures for all cases
    
    if case_list is None:
        # Default to run all cases
        case_list = ["Flat","A","B","C","D","E","F"]
    
    # Convert case_list to list if necessary
    if type(case_list) is not list:
        case_list = [case_list]
    
    # Define members to investigate for each case
    titles = dict()
    titles['Flat']       = ["Primary Members", "Secondary Members", "Total Load"]
    titles['A']          = ["Top Primary Member", "Secondary Members"]
    titles['B']          = ["Top Primary Member", "Bottom Primary Member", "Secondary Members"]
    titles['C']          = ["Primary Members", "Secondary Member 2"]
    titles['D']          = ["Primary Members", "Secondary Member 1"]
    titles['E']          = ["Top Primary Member", "Bottom Primary Member", "Secondary Member 1"]
    titles['F']          = ["Top Primary Member", "Secondary Member 2"]
    
    Cs_list = [0.001,0.1,0.2,0.3]

    num_spaces_list = [2,16]
            
    # Loop over cases, Cs, num_spaced, and member_title to make figures        
    for case in case_list:
        for Cs in Cs_list:
            for num_spaces in num_spaces_list:
                for member_title in titles[case]:
                    two_way_plot(case, Cs, num_spaces, member_title)

if __name__ == "__main__":

    try:
        os.mkdir(os.path.join('figures'))
    except FileExistsError:
        pass  

    create_paper_figures()
    create_all_figures()
    #two_way_plot("A",0.3,16,"Top Primary Member","Test_Plot")
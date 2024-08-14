from PyPonding.structures import IdealizedBay
from math import pi, cos, cosh
import numpy as np

# Based on Example 1 from Ruddy, J. L. (1986). “Ponding of Concrete Deck Floors.” Engineering Journal, AISC, 23(3), 107–115.

# Define units
inch    = 1.0
kip     = 1.0
lb      = kip/1000.0
ft      = 12.0*inch
in_per_ft = inch/ft
ksi     = kip/inch**2
plf     = lb/ft
psf     = lb/ft**2
pcf     = psf/ft
kipft   = kip*ft
cy      = (36*inch)**3

# Bay Geometry
Lp = 28*ft
Ls = 28*ft
num_spaces = 4
S  = Lp/num_spaces

# Loads
qD = 49*psf
γw = 145*pcf   
zw = 0 # Concrete above undeflected floor treated as dead load

# Top of roof elevation
zTL = 0.0*inch
zTR = 0.0*inch
zBL = 0.0*inch
zBR = 0.0*inch

# Material properties
E = 29000.0*ksi

# Member stiffness
Ip = 984*inch**4
Is = 199*inch**4

# Define IdealizedBay object
input = {
    'primary_member_span': Lp,
    'secondary_member_span': Ls,
    'number_of_joist_spaces': num_spaces,
    'dead_load_uniform': qD,
    'dead_load_primary_member': 0*plf,
    'water_density': γw,  
    'snow_density': 0.0*pcf,
    'snow_height': 0.0*inch,
    'alpha': 1.0,
    'load_factor_dead':    1.0,
    'load_factor_ponding': 1.0,
    'load_factor_snow':    0.0,
    'z_TL': zTL,
    'z_TR': zTR,
    'z_BL': zBL,
    'z_BR': zBR,
    'secondary_member_camber': 0.000*inch,
    'primary_member_camber_T': 0.000*inch,
    'primary_member_camber_B': 0.000*inch,
    'edge_condition_L': 'mirrored',
    'edge_condition_R': 'mirrored',
    'edge_condition_T': 'mirrored',
    'edge_condition_B': 'mirrored',
    'E':  E,
    'As': 100*inch**2,
    'Ap': 100*inch**2,
    'Is': Is,
    'Ip': Ip,
    'analsis_engine': 'FE',
}
bay = IdealizedBay(**input)


bay.include_ponding_effect = False
results0 = bay.Run_Analysis(zw)
bay.include_ponding_effect = True
resultsP = bay.Run_Analysis(zw)

print('==== Total Load ====')
print(f'No Ponding Effect:    {results0.total_factored_load:.3f} kips')
print(f'With Ponding Effect:  {resultsP.total_factored_load:.3f} kips')
print(f'Amplification Factor: {resultsP.total_factored_load/results0.total_factored_load:.3f}')

print('==== Additional Volume of Concrete ====')
additional_volume = (resultsP.total_factored_load-qD*Lp*Ls)/γw
print(f'{additional_volume/cy:.3f} yd^3')

from math import pi
from sji_load_tables import get_joist_data, lightest_joist
from libdenavit import JoistGirder


# Units
inch = 1.0
kip = 1.0
lb  = kip/1000.0
ft  = 12.0*inch
in_per_ft = inch/ft
ksi = kip/inch**2
plf = lb/ft
psf = lb/ft**2
pcf = psf/ft
kipft = kip*ft

# Bay Geometry
Lp = 50*ft
Ls = 60*ft
S = Lp/8
E = 29000*ksi

# Loads
D  = 15*psf
Lr = 20*psf

w_D = D*S
w_Lr = Lr*S
w_R = 169*plf # Figure 18

P_D = D*S*Ls
P_Lr = Lr*S*Ls
P_R = 1444*lb

gamma_w = 62.4*pcf

# Method Set 1 - No Rain
print('\nMethod Set 1')
joist = lightest_joist(Ls/ft,(w_D+w_Lr)/plf,0,design_basis='ASD')
print((w_D+w_Lr)/plf)
print(joist)
print(f'{P_D+P_Lr} kips')

# Method Set 2 - Rain, No Ponding
print('\nMethod Set 2')
joist = lightest_joist(Ls/ft,(w_D+w_R)/plf,0,design_basis='ASD')
print(w_D/plf)
print(w_R/plf)
print((w_D+w_R)/plf)
print(joist)
P_TL = max(P_D+P_Lr,P_D+P_R)
print(f'{P_TL} kips')

# Method Set 3 - Rain, Ponding using 2019 AF
print('\nMethod Set 3')
joist_designation = '40LH10' #'44LH11'
P_TL = 13.2

joist = get_joist_data(joist_designation,Ls/ft)
print(joist)
Is = joist.approx_moment_of_inertia(shear_deformation_factor=1.15)

joist_girder = JoistGirder('ASD',Lp/ft,48,8,P_TL)
Ip = joist_girder.moment_of_inertia()/1.15

Cs = (gamma_w*S*Ls**4)/(pi**4*E*Is)
Cp = (gamma_w*Ls*Lp**4)/(pi**4*E*Ip)
Bp = 1/(1-1.15*Cp-Cs)

print(f'Is = {Is:.3f}')
print(f'Ip = {Ip:.3f}')
print(f'Cs = {Cs:.3f}')
print(f'Cp = {Cp:.3f}')
print(f'Bp = {Bp:.3f}')

P_TL_calc = max(P_D+P_Lr,Bp*(P_D+P_R))
if P_TL >= P_TL_calc:
    print(f'{P_TL:.2f} >= {P_TL_calc:.2f} --- OKAY')
else:
    print(f'{P_TL:.2f} < {P_TL_calc:.2f} --- NO GOOD')

total_load = joist.total_load('ASD')
total_load_calc = max(w_D+w_Lr,Bp*(w_D+w_R))/plf
if total_load >= total_load_calc:
    print(f'{total_load:.2f} >= {total_load_calc:.2f} --- OKAY')
else:
    print(f'{total_load:.2f} < {total_load_calc:.2f} --- NO GOOD')
    
    
# Method Set 4 - Rain, Ponding using Proposed AF
print('\nMethod Set 4')
joist_designation = '44LH09'
P_TL = 13.2

joist = get_joist_data(joist_designation,Ls/ft)
print(joist)
Is = joist.approx_moment_of_inertia(shear_deformation_factor=1.15)

joist_girder = JoistGirder('ASD',Lp/ft,48,8,P_TL)
Ip = joist_girder.moment_of_inertia()/1.15

Cs = (gamma_w*S*Ls**4)/(pi**4*E*Is)
Cp = (gamma_w*Ls*Lp**4)/(pi**4*E*Ip)

print(f'Cs = {Cs:.3f}')
print(f'Cp = {Cp:.3f}')

zw_over_zh = 10/15
Bpmaxp = 1/(1-0.60*Cp-1.15*Cs)
Bpmaxs = 1/(1-0.55*Cp-1.00*Cs)
if zw_over_zh <= 0.20:
    Bpp = 1.00
    Bps = 1.00
elif zw_over_zh <= 0.85:
    Bpp = 1.00 + (Bpmaxp - 1.00)*(zw_over_zh-0.20)/(0.85-0.20)
    Bps = 1.00 + (Bpmaxs - 1.00)*(zw_over_zh-0.20)/(0.85-0.20)
else:
    Bpp = Bpmaxp
    Bps = Bpmaxs

print(f'Bps = {Bps:.3f}')
print(f'Bpp = {Bpp:.3f}')

P_TL_calc = max(P_D+P_Lr,Bpp*(P_D+P_R))
if P_TL >= P_TL_calc:
    print(f'{P_TL:.2f} >= {P_TL_calc:.2f} --- OKAY')
else:
    print(f'{P_TL:.2f} < {P_TL_calc:.2f} --- NO GOOD')

total_load = joist.total_load('ASD')
total_load_calc = max(w_D+w_Lr,Bps*(w_D+w_R))/plf
if total_load >= total_load_calc:
    print(f'{total_load:.2f} >= {total_load_calc:.2f} --- OKAY')
else:
    print(f'{total_load:.2f} < {total_load_calc:.2f} --- NO GOOD')
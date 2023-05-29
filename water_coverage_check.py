import numpy as np
import matplotlib.pyplot as plt

# Creates figures that show the coverage of water on a bay for Case E and Case F
# Assumes no deflections and no camber

num_secondary_members = 7
zw_over_zh = 0.5


# CASE E
secondary_member_y_positions = np.linspace(0,1,num_secondary_members).tolist()
secondary_member_z_positions = np.linspace(0,1,num_secondary_members).tolist()
water_y_positions = np.linspace(0,1,num_secondary_members).tolist()

water_x_positions = []
for i in range(num_secondary_members):
    if secondary_member_z_positions[i] <= zw_over_zh:
        x = (zw_over_zh-secondary_member_z_positions[i])/(1-secondary_member_z_positions[i])
        water_x_positions.append(x)
    else:
        water_x_positions.append(0)

plt.figure(figsize=(3.25,3.25))        
plt.fill_between(water_x_positions,water_y_positions)
for member in range(num_secondary_members):
    plt.plot([0,1],[secondary_member_y_positions[member],secondary_member_y_positions[member]],'k')
plt.ylim(0,1)
plt.xlim(0,1)
plt.xticks([])
plt.yticks([])
plt.title('Case E')


# CASE F
secondary_member_y_positions = np.linspace(0,1,num_secondary_members).tolist()
secondary_member_z_positions = np.linspace(1,0,num_secondary_members).tolist()
water_y_positions = np.linspace(0,1,num_secondary_members).tolist()

water_x_positions = []
for i in range(num_secondary_members):
    if secondary_member_z_positions[i] <= zw_over_zh:
        x = 1- (zw_over_zh-secondary_member_z_positions[i])/(1-secondary_member_z_positions[i])
        water_x_positions.append(x)
    else:
        water_x_positions.append(1)

water_x_positions.append(0)
water_y_positions.append(1)

plt.figure(figsize=(3.25,3.25))        
plt.fill_between(water_x_positions,water_y_positions)
for member in range(num_secondary_members):
    plt.plot([0,1],[secondary_member_y_positions[member],secondary_member_y_positions[member]],'k')
plt.ylim(0,1)
plt.xlim(0,1)
plt.xticks([])
plt.yticks([])
plt.title('Case F')


plt.show()
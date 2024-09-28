# ecef_to_sez.py
#
# Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
# Converts SEZ frame parameters with respect to the ECEF origin to SEZ components
# which consists of a translation and two inverse rotations from ECEF vector
# Parameters:
#  o_x_km: x-component of SEZ with respect to ECEF origin
#  o_y_km: y-component of SEZ with respect to ECEF origin
#  o_z_km: z-component of SEZ with respect to ECEF origin
#  x_km: ECEF x-position 
#  y_km: ECEF y-position 
#  z_km: ECEF z-position 
# Output:
#  s_km: South position in SEZ
#  e_km: East position in SEZ
#  z_km: Height position in SEZ
#
# Written by Mandar Ajgaonkar
# Other contributors: None
#
# Optional license statement, e.g., See the LICENSE file for the license.

# import Python modules
import math  # math module
import sys   # argv

# "constants"
R_E_KM = 6378.1363 # radius of Earth in km
E_E = 0.081819221456 # unitless

# helper functions

# function description

## calculate denominator 
def calc_denom(ecc, lat_rad):
    return math.sqrt(1.0 - (ecc**2.0) * (math.sin(lat_rad))**2)

# initialize script arguments
o_x_km = float('nan')  
o_y_km = float('nan')  
o_z_km = float('nan')  
x_km = float('nan')
y_km = float('nan')
z_km = float('nan')

# parse script arguments 
if len(sys.argv) == 7:
    try:
        o_x_km = float(sys.argv[1])
        o_y_km = float(sys.argv[2])
        o_z_km = float(sys.argv[3])
        x_km = float(sys.argv[4])
        y_km = float(sys.argv[5])
        z_km = float(sys.argv[6])
    except ValueError:
        print("Error: o_x_km, o_y_km, o_z_km, x_km, y_km, and z_km must be numeric.")
        exit()
else:
    print('python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km')
    exit()

# calculate ECEF vector from the station to object / satellite
x_ECEF_km = x_km - o_x_km
y_ECEF_km = y_km - o_y_km
z_ECEF_km = z_km - o_z_km

# calculate longitude
lon_rad = math.atan2(o_y_km,o_x_km)
lon_deg = lon_rad*180.0/math.pi

# initialize lat_rad, r_lon_km, r_z_km
lat_rad = math.asin(o_z_km/math.sqrt(o_x_km**2+o_y_km**2+o_z_km**2))
r_lon_km = math.sqrt(o_x_km**2+o_y_km**2)
prev_lat_rad = float('nan')

# iteratively find latitude
c_E = float('nan')
count = 0
while (math.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad)>10e-7) and count<5:
  denom = calc_denom(E_E,lat_rad)
  c_E = R_E_KM/denom
  prev_lat_rad = lat_rad
  lat_rad = math.atan((o_z_km+c_E*(E_E**2)*math.sin(lat_rad))/r_lon_km)
  count = count+1
  
# calculate hae
hae_km = r_lon_km/math.cos(lat_rad)-c_E

# multiply r_SEZ by Rz^-1(lon_rad)
vec_x = x_ECEF_km*math.cos(lon_rad) + y_ECEF_km*math.sin(lon_rad)
vec_y = x_ECEF_km*-math.sin(lon_rad) + y_ECEF_km*math.cos(lon_rad)
vec_z = z_ECEF_km 

# multiply above vector by Ry^-1(90 - lat_rad) to get s_km, e_km and z_km 
s_km = vec_x*math.sin(lat_rad) - vec_z*math.cos(lat_rad)
e_km = vec_y
z_km = vec_x*math.cos(lat_rad) + vec_z*math.sin(lat_rad)

# print results
print(s_km)
print(e_km)
print(z_km)

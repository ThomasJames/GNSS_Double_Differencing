from math import sqrt, cos, sin, degrees, acos
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from numpy import transpose, linalg

"""
The data contains 1 of the three epochs of phase and pseudorage observations measured on a calibration baseline in valencia, spain.
The sensors used are geodetic quality recievers using choke ring antennas.
The reciever on pillar 1A is treated as the reference reciever.
The reciever on pillar 3A is treated as the monument.
"""

# X, Y, Z ECEF coordinates for the phase center of the receiver
pillar_1A_base = np.array([[4929635.400], [-29041.877], [4033567.846]])  # Reference receiver
pillar_3A_rover = np.array([[4929605.400], [-29123.700], [4033603.800]])  # Monument

distance_between_receivers = 94.4  # Measured in meters / Approximate

"""
ECEF coordinates (m) of the satelite phase centers when they transmitted the signals measured at each epoch.
"""
# ECEF SATELLITE POSITIONS X, Y, Z (m) (already corrected for earth rotation during signal travel time)
# 2016 11 15 22 19  5

G10 = [4634093.207, -19899701.050, 16933747.321]
G12 = [22559170.178, -8979632.676, 10377257.530]
G13 = [23277536.897, 12575815.276, -2029027.200]
G15 = [25950462.808, 2443858.353, 5881092.070]
G17 = [5785091.956, 16827408.400, 20125597.869]
G18 = [13564948.214, -21357948.777, 8232124.013]
G19 = [12262838.101, 17165601.305, 15682863.092]
G24 = [15569244.807, -1039249.482, 21443791.252]

"""
Use double differenced phase measurements, from the first epoch of data only 2016_11_15_22_19_5
TO compute the precise coordinates of the pillar 3A sensor phase center. 
These are pseudo range measurements 
"""

# BASE OBSERVATIONS (Pillar 1A) C1C (metres), L1C (L1 cycles)
# 2016_11_15_22_19_5

G24_base_obs = [20436699.926, 107395596.426]
G19_base_obs = [22181729.713, 116565751.296]
G18_base_obs = [23217805.737, 122010370.583]
G17_base_obs = [23379753.757, 122861442.012]
G15_base_obs = [21346539.664, 112176830.803]
G13_base_obs = [23087780.798, 121327099.499]
G12_base_obs = [20647534.024, 108503516.027]
G10_base_obs = [23726969.123, 124686036.295]

# ROVER OBSERVATIONS (Pillar 3A) C1C (metres), L1C (L1 cycles)
# 2016 11 15 22 19  5 (First Epoch)
# meters, cycles

G24_rover_obs = [20436682.002, 107395502.123]
G19_rover_obs = [22181785.598, 116566080.299]
G18_rover_obs = [23217736.821, 122010019.631]
G17_rover_obs = [23379790.820, 122861635.973]
G15_rover_obs = [21346576.786, 112177022.660]
G13_rover_obs = [23087860.345, 121327512.345]
G12_rover_obs = [20647514.655, 108503447.644]
G10_rover_obs = [23726881.094, 124685588.685]

# At the first epoch we have 16 raw phase observations in cycles.

l = np.transpose([
    G24_base_obs[1], G24_rover_obs[1],
    G19_base_obs[1], G19_rover_obs[1],
    G18_base_obs[1], G18_rover_obs[1],
    G17_base_obs[1], G17_rover_obs[1],
    G15_base_obs[1], G15_rover_obs[1],
    G13_base_obs[1], G13_rover_obs[1],
    G12_base_obs[1], G12_rover_obs[1],
    G10_base_obs[1], G10_rover_obs[1]
])


"""
Phase Ambiguity terms (N) for each measurement, before and after ambiguity resolution 
Use integer terms in computations
2016 11 15 22 19  5
G24 is a reference satelite - and has the highest 
"""

# Phase ambiguities for each epoch and each phase measurement:

before_ambiguity_resolution = np.array([[4929605.364], [-29123.817], [4033603.867]])
G24toG10_before = 12.564
G24toG12_before = 34.873
G24toG13_before = -3.838
G24toG15_before = -4.170
G24toG17_before = 1.538
G24toG18_before = 11.324
G24toG19_before = 34.352

after_ambiguity_resolution = np.array([[4929605.542], [-29123.828], [4033603.932]])
G24toG10_after = int(12.000)
G24toG12_after = int(35.000)
G24toG13_after = int(-4.000)
G24toG15_after = int(-4.000)
G24toG17_after = int(1.000)
G24toG18_after = int(11.000)
G24toG19_after = int(34.000)

G24toG10_noise = G24toG10_before - G24toG10_after
G24toG12_noise = G24toG12_before - G24toG12_after
G24toG13_noise = G24toG13_before - G24toG13_after
G24toG15_noise = G24toG15_before - G24toG15_after
G24toG17_noise = G24toG17_before - G24toG17_after
G24toG18_noise = G24toG18_before - G24toG18_after
G24toG19_noise = G24toG19_before - G24toG19_after

# 16 x 8:  Differencing matrix
S = np.array([[1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1]

              ])

D = np.array([[1, -1, 0, 0, 0, 0, 0, 0],
              [1, 0, -1, 0, 0, 0, 0, 0],
              [1, 0, 0, -1, 0, 0, 0, 0],
              [1, 0, 0, 0, -1, 0, 0, 0],
              [1, 0, 0, 0, 0, -1, 0, 0],
              [1, 0, 0, 0, 0, 0, -1, 0],
              [1, 0, 0, 0, 0, 0, 0, -1]])

"""
How to compute the signal wavelength:
λ=𝑐/𝑓
c = speed of light in vacuum (299792458.0 ms-1)
f = signal frequency (L1: 1575.42MHz, L2: 1227.6MHz)
"""

c = 299792458.0
f = 1575.42

wl = c / f

"""
Standard Deviations
L1C variance for satellite at elevation angle E: s*s/sin(E)
NOTE: Satellite G24 has the highest elevation: 71 degrees.

"""

# Calculate Satelite elevation from a local horizon.
def elevation_calculator(ECtoSat_Z, ECtoSat_X, EctoRec_Z, EctoRec_X, SattoRec_range):

    EC_sat_range = sqrt((ECtoSat_X ** 2) + (ECtoSat_Z ** 2))
    EctoRec_range = sqrt((EctoRec_X) ** 2 + (EctoRec_Z ** 2))
    SattoRec_range = SattoRec_range

    c = EC_sat_range
    b = EctoRec_range
    a = SattoRec_range


    angle = degrees((acos((c ** 2 - b ** 2 - a ** 2) / (-2.0 * a * b))))

    elevation = abs(angle - 90)
    return elevation


def variance(s, e):
    a = s**2 * s / cos(e)
    return a

l1_SD = 0.003

# Determine elevationss
G24_elevation = (elevation_calculator(G24[2], G24[0], pillar_1A_base[2], pillar_1A_base[0], G24_base_obs[0]))
G19_elevation = (elevation_calculator(G19[2], G19[0], pillar_1A_base[2], pillar_1A_base[0], G19_base_obs[0]))
G17_elevation = (elevation_calculator(G17[2], G17[0], pillar_1A_base[2], pillar_1A_base[0], G17_base_obs[0]))
G15_elevation = (elevation_calculator(G15[2], G15[0], pillar_1A_base[2], pillar_1A_base[0], G15_base_obs[0]))
G13_elevation = (elevation_calculator(G13[2], G13[0], pillar_1A_base[2], pillar_1A_base[0], G13_base_obs[0]))
G12_elevation = (elevation_calculator(G12[2], G12[0], pillar_1A_base[2], pillar_1A_base[0], G12_base_obs[0]))
G10_elevation = (elevation_calculator(G10[2], G10[0], pillar_1A_base[2], pillar_1A_base[0], G10_base_obs[0]))
G18_elevation = 5 # Something wrong with this satellite..(Random Value for now)

# Populate the a vector of variances.

G24_variance = variance(l1_SD, G24_elevation)
G19_variance = variance(l1_SD, G19_elevation)
G17_variance = variance(l1_SD, G17_elevation)
G15_variance = variance(l1_SD, G15_elevation)
G13_variance = variance(l1_SD, G13_elevation)
G12_variance = variance(l1_SD, G12_elevation)
G10_variance = variance(l1_SD, G10_elevation)
G18_variance = variance(l1_SD, G18_elevation)

# Vector of variances
variances = np.array([
[G24_variance],
[G19_variance],
[G17_variance],
[G15_variance],
[G13_variance],
[G12_variance],
[G10_variance],
[G18_variance],
])




# print("G24 angle of elevation:", G24_elevation)
# print("G19 angle of elevation:", G19_elevation)
# print("G17 angle of elevation:", G17_elevation)
# print("G15 angle of elevation:", G15_elevation)
# print("G13 angle of elevation:", G13_elevation)
# print("G12 angle of elevation:", G12_elevation)
# print("G10 angle of elevation:", G10_elevation)
# print("G18 angle of elevation:")














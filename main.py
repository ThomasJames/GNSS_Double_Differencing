from math import sqrt, cos, sin, degrees, acos
import numpy as np
from numpy import transpose, linalg
from Computations import DD, Variance, matrix_heatmap, vector_heatmap, flipped_vector_heatmap, MatrixOperations
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List

"""
GOAL: Calculate the coordinates of the reference antenna (ARP) of the roving receiver 

Try to get your answer close to the figures for 3A. The nominal coordinates given mean you do not need to iterate the 
least squares solution, you should converge on the answer with on round of matrix inversion
The data contains 1 of the three epochs of phase and pseudorange observations measured on a calibration baseline in 
valencia, spain.
The sensors used are geodetic quality receivers using choke ring antennas.
Pillar 1A is treated as the reference receiver.
Pillar 3A is treated as the monument.
"""

# X, Y, Z ECEF coordinates for the phase center of the receiver
pillar_1A_base = np.array([[4929635.440], [-29041.877], [4033567.846]])  # Reference receiver

# Trying to reproduce these coordinates
pillar_3A_rover = np.array([[4929605.400], [-29123.700], [4033603.800]])  # Monument
distance_between_receivers = 94.4  # Measured in meters / Approximate
l1_SD = 0.003

"""
ECEF coordinates (m) of the satellite phase centers when they transmitted the signals measured at each epoch.
"""

# Make a note of the L1 Wl
wl = 0.19029367

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
    G24_base_obs[1],
    G24_rover_obs[1],
    G19_base_obs[1],
    G19_rover_obs[1],
    G18_base_obs[1],
    G18_rover_obs[1],
    G17_base_obs[1],
    G17_rover_obs[1],
    G15_base_obs[1],
    G15_rover_obs[1],
    G13_base_obs[1],
    G13_rover_obs[1],
    G12_base_obs[1],
    G12_rover_obs[1],
    G10_base_obs[1],
    G10_rover_obs[1]])

"""
Phase Ambiguity terms (N) for each measurement, before and after ambiguity resolution 
Use integer terms in computations
2016 11 15 22 19  5
G24 is a reference satellite - and has the highest 
"""

# Phase ambiguities for each epoch and each phase measurement:
before_ambiguity_resolution = np.array([[4929605.364], [-29123.817], [4033603.867]])
G24toG19_before = 34.352
G24toG18_before = 11.324
G24toG17_before = 1.538
G24toG15_before = -4.170
G24toG13_before = -3.838
G24toG12_before = 34.873
G24toG10_before = 12.564

after_ambiguity_resolution = np.array([[4929605.542], [-29123.828], [4033603.932]])
G24toG19_after = 34.000
G24toG18_after = 11.000
G24toG17_after = 1.000
G24toG15_after = -4.000
G24toG13_after = -4.000
G24toG12_after = 35.000
G24toG10_after = 12.000

a_a_r = [G24toG19_after,
         G24toG18_after,
         G24toG17_after,
         G24toG15_after,
         G24toG13_after,
         G24toG12_after,
         G24toG10_after]

G24_base_var = Variance(sat_coords=G24, receiver_coords=pillar_1A_base, L1=True)
G24_rover_var = Variance(sat_coords=G24, receiver_coords=pillar_3A_rover, L1=True)
G19_base_var = Variance(sat_coords=G19, receiver_coords=pillar_1A_base, L1=True)
G19_rover_var = Variance(sat_coords=G19, receiver_coords=pillar_3A_rover, L1=True)
G18_base_var = Variance(sat_coords=G18, receiver_coords=pillar_1A_base, L1=True)
G18_rover_var = Variance(sat_coords=G18, receiver_coords=pillar_3A_rover, L1=True)
G17_base_var = Variance(sat_coords=G17, receiver_coords=pillar_1A_base, L1=True)
G17_rover_var = Variance(sat_coords=G17, receiver_coords=pillar_3A_rover, L1=True)
G15_base_var = Variance(sat_coords=G15, receiver_coords=pillar_1A_base, L1=True)
G15_rover_var = Variance(sat_coords=G15, receiver_coords=pillar_3A_rover, L1=True)
G13_base_var = Variance(sat_coords=G13, receiver_coords=pillar_1A_base, L1=True)
G13_rover_var = Variance(sat_coords=G13, receiver_coords=pillar_3A_rover, L1=True)
G12_base_var = Variance(sat_coords=G12, receiver_coords=pillar_1A_base, L1=True)
G12_rover_var = Variance(sat_coords=G12, receiver_coords=pillar_3A_rover, L1=True)
G10_base_var = Variance(sat_coords=G10, receiver_coords=pillar_1A_base, L1=True)
G10_rover_var = Variance(sat_coords=G10, receiver_coords=pillar_3A_rover, L1=True)

elevations_radians = np.array([
    [G24_base_var.elevation_calculator()],
    [G24_rover_var.elevation_calculator()],
    [G19_base_var.elevation_calculator()],
    [G19_rover_var.elevation_calculator()],
    [G18_base_var.elevation_calculator()],
    [G18_rover_var.elevation_calculator()],
    [G17_base_var.elevation_calculator()],
    [G17_rover_var.elevation_calculator()],
    [G15_base_var.elevation_calculator()],
    [G15_rover_var.elevation_calculator()],
    [G13_base_var.elevation_calculator()],
    [G13_rover_var.elevation_calculator()],
    [G12_base_var.elevation_calculator()],
    [G12_rover_var.elevation_calculator()],
    [G10_base_var.elevation_calculator()],
    [G10_rover_var.elevation_calculator()]])

satelltie_names = np.array([["Base to G19"],
                            ["Rover to G19"],
                            ["Base to  G18"],
                            ["Rover to G18"],
                            ["Base to G17"],
                            ["Rover to G17"],
                            ["Base to G15"],
                            ["Rover to G15"],
                            ["Base to G13"],
                            ["Rover to G13"],
                            ["Base to G12"],
                            ["Rover to G12"],
                            ["Base to G10"],
                            ["Rover to G10"],
                            ["Base to G24"],
                            ["Rover to G24"]])

# Elevations in degrees
elevations_degrees = np.array([degrees(x) for x in elevations_radians])

variance_vector = np.array([
    [G24_base_var.variance()],
    [G24_rover_var.variance()],
    [G19_base_var.variance()],
    [G19_rover_var.variance()],
    [G18_base_var.variance()],
    [G18_rover_var.variance()],
    [G17_base_var.variance()],
    [G17_rover_var.variance()],
    [G15_base_var.variance()],
    [G15_rover_var.variance()],
    [G13_base_var.variance()],
    [G13_rover_var.variance()],
    [G12_base_var.variance()],
    [G12_rover_var.variance()],
    [G10_base_var.variance()],
    [G10_rover_var.variance()]])

# 16 x 8:  Differencing matrix
S = np.array([[1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1]])

D = np.array([[1, -1, 0, 0, 0, 0, 0, 0],
              [1, 0, -1, 0, 0, 0, 0, 0],
              [1, 0, 0, -1, 0, 0, 0, 0],
              [1, 0, 0, 0, -1, 0, 0, 0],
              [1, 0, 0, 0, 0, -1, 0, 0],
              [1, 0, 0, 0, 0, 0, -1, 0],
              [1, 0, 0, 0, 0, 0, 0, -1]])


if __name__ == "__main__":

    flipped_vector_heatmap(variance_vector, "Variances")
    flipped_vector_heatmap(elevations_radians, "Radians")
    flipped_vector_heatmap(elevations_degrees, "Degrees")

    """
    l vector - The vector of raw observations
    """

    flipped_vector_heatmap(l, "Vector of observations (l)")

    """
    For purposes of single differencing and double differencing.
    The D and S are created.
    D DIM: 7 x 8
    S DIM: 8 X 16
    """

    matrix_heatmap(S, "Single Differencing (S)")

    matrix_heatmap(D, "Double Differencing (D)")

    """
    SINGLE DIFFERENCING
    As receiver 1A and 3A are relatively close together (9.4m)
    The ionspheric and tropospheric terms cancel out
    Therefore we can calculate receiver-receiver single difference (RRSD)
    sl = vector of receiver-receiver single differences 
    DIM: 8 x 1
    """

    sl = S.dot(l)

    flipped_vector_heatmap(sl, "Vector of single differences (sl)")

    """
    DOUBLE DIFFERENCING 
    
    We choose the satellite with the highest elevation as the reference satellite.
    This will be G24 - This satellite will have the least noise and multipath 
    
    Dsl = vector of double differences
    DIM: 7 x 1
    """

    Dsl = D.dot(sl)

    """
    Subtract the corresponding ambiguity resolved phase ambiguity term. 
    """
    DD_s_p_a = []
    for i in range(len(Dsl)):
        DD_s_p_a.append(Dsl[i] - a_a_r[i])

    flipped_vector_heatmap(Dsl, "Double differences (Dsl)")

    """
    Covariance matrix of the observation vector.
    This is an identity matrix scaled by the variance.
    It is assumed that the variance of each raw phase observation has the same l1 standard deviation is 0.003. 
    Therefore variance is this value squared.
    UNITS: Cycles
    DIM: 16 x 16 
    """

    cl = np.eye(16, 16) * (1 / wl * variance_vector)  # back to cycles
    matrix_heatmap(cl, "Covariance matrix of observations (cl)")

    """
    Obtain the covariance matrix of the double differences.
    This is because we require Wd - The weight matrix of the double difference vector. (The inverse of cd)
    Apply gauss's propagation of error law: y = Qx, Cx to work out Cd
    Where Q is the matrix operator - No uncertainties associated
    Where x is stochastic quantity
    Where y is stochastic
    Where Cx is the covariance matrix of x and is known
    
    d = DSl 
    Cl Exists and is known.
    
    Use the formula: Cd = DSCl (DS)^T 
    
    DIM: 7 x 7 
    """

    Cd = MatrixOperations(D=D, S=S, Cl=cl)
    Cd = Cd.Cd_calculator()
    matrix_heatmap(Cd, "covariance matrix of the double differences (Cd)")

    """
    Calculate the weight matrix of the double differences. 
    
    Wd = Cd^-1
    DIM: 7 x 7
    """

    Wd = linalg.inv(Cd)
    matrix_heatmap(Wd, "Weight (Wd)")

    """

    """
    G24G19 = DD(ref_station=pillar_1A_base, rov_station=pillar_3A_rover, corresponding_sat=G19, sat_ref=G24,
                observed=DD_s_p_a[0])

    G24G18 = DD(ref_station=pillar_1A_base, rov_station=pillar_3A_rover, corresponding_sat=G18, sat_ref=G24,
                observed=DD_s_p_a[1])

    G24G17 = DD(ref_station=pillar_1A_base, rov_station=pillar_3A_rover, corresponding_sat=G17, sat_ref=G24,
                observed=DD_s_p_a[2])

    G24G15 = DD(ref_station=pillar_1A_base, rov_station=pillar_3A_rover, corresponding_sat=G15, sat_ref=G24,
                observed=DD_s_p_a[3])

    G24G13 = DD(ref_station=pillar_1A_base, rov_station=pillar_3A_rover, corresponding_sat=G13, sat_ref=G24,
                observed=DD_s_p_a[4])

    G24G12 = DD(ref_station=pillar_1A_base, rov_station=pillar_3A_rover, corresponding_sat=G12, sat_ref=G24,
                observed=DD_s_p_a[5])

    G24G10 = DD(ref_station=pillar_1A_base, rov_station=pillar_3A_rover, corresponding_sat=G10, sat_ref=G24,
                observed=DD_s_p_a[6])

    """
    Calculate the b vector:
    This is the observed double differencing measurements - the computed.  
    """
    b = np.array([[G24G19.calc_b_vector()],
                  [G24G18.calc_b_vector()],
                  [G24G17.calc_b_vector()],
                  [G24G15.calc_b_vector()],
                  [G24G13.calc_b_vector()],
                  [G24G12.calc_b_vector()],
                  [G24G10.calc_b_vector()]])

    vector_heatmap(b, "Observed - Computed")

    # Populate the design matrix
    A = np.array([[G24G19.x_diff(), G24G19.y_diff(), G24G19.z_diff()],
                  [G24G18.x_diff(), G24G18.y_diff(), G24G18.z_diff()],
                  [G24G17.x_diff(), G24G17.y_diff(), G24G17.z_diff()],
                  [G24G15.x_diff(), G24G15.y_diff(), G24G15.z_diff()],
                  [G24G13.x_diff(), G24G13.y_diff(), G24G13.z_diff()],
                  [G24G12.x_diff(), G24G12.y_diff(), G24G12.z_diff()],
                  [G24G10.x_diff(), G24G10.y_diff(), G24G10.z_diff()]])

    matrix_heatmap(A, "Design (A)")

    """
    Output the ATWA matrix 
    """
    atwa = MatrixOperations(A=A, W=Wd)
    atwa = atwa.ATWA()
    matrix_heatmap(atwa, "ATWA")

    """
    Output the (ATWA)^-1 matrix 
    """
    inverse_atwa = linalg.inv(atwa)

    matrix_heatmap(inverse_atwa, "(ATWA)^-1")

    ATWb = (transpose(A).dot(Wd)).dot(b)

    x_hat = inverse_atwa.dot(ATWb)

    x, y, z = x_hat
    print(x)
    print(y)
    print(z)

    print(f"Updated Cooridnates: {pillar_3A_rover[0] + x} actual coordinates: {after_ambiguity_resolution[0]}")
    print(f"Updated Cooridnates: {pillar_3A_rover[1] + y} actual coordinates: {after_ambiguity_resolution[1]}")
    print(f"Updated Cooridnates: {pillar_3A_rover[2] + z} actual coordinates: {after_ambiguity_resolution[2]}")


    # Quality Assessment
    # Distance between nominal reciever and reference receiver
    def distance(point_1: List[float], point_2: List[float]) -> float:
        """""
        Find the difference between two points given sets of [X, Y, Z] coordinates.
        """""
        return sqrt((point_2[0] - point_1[0]) ** 2 +
                    (point_2[1] - point_1[1]) ** 2 +
                    (point_2[2] - point_1[2]) ** 2)


    print(distance(pillar_1A_base, pillar_3A_rover))
    print(distance(pillar_1A_base, pillar_3A_rover + x_hat))

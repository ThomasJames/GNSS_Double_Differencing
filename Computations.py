from math import sqrt, cos, sin, degrees, acos
import numpy as np
from numpy import transpose, linalg
from Matrix_Computation_Classes import DD, Variance, distance
import matplotlib.pyplot as plt
from matplotlib import transforms
import seaborn as sns
import pandas as pd
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
pillar_1A_base = np.array([[4929635.400], [-29041.877], [4033567.846]])  # Reference receiver

# Trying to reproduce these coordinates
pillar_3A_rover = np.array([[4929605.400], [-29123.700], [4033603.800]])  # Monument
distance_between_receivers = 94.4  # Measured in meters / Approximate

l1_SD = 0.003

"""
ECEF coordinates (m) of the satellite phase centers when they transmitted the signals measured at each epoch.
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

# Calculate the noise of each measurement.
G24toG19_noise = G24toG19_before - G24toG19_after
G24toG18_noise = G24toG18_before - G24toG18_after
G24toG17_noise = G24toG17_before - G24toG17_after
G24toG15_noise = G24toG15_before - G24toG15_after
G24toG13_noise = G24toG13_before - G24toG13_after
G24toG12_noise = G24toG12_before - G24toG12_after
G24toG10_noise = G24toG10_before - G24toG10_after


G24_base_var = Variance(sat_coords=G24, receiver_coords=pillar_1A_base, L1=True)
G24_rover_var = Variance(sat_coords=G24, receiver_coords=pillar_3A_rover, L1=True)
G19_base_var = Variance(sat_coords=G19, receiver_coords=pillar_1A_base, L1=True)
G19_rover_var = Variance(sat_coords=G19, receiver_coords=pillar_3A_rover, L1=True)
G18_base_var = Variance(sat_coords=G18, receiver_coords=pillar_1A_base, L1=True)
G18_rover_var = Variance(sat_coords=G18, receiver_coords=pillar_3A_rover, L1=True)
G17_base_var = Variance(sat_coords=G17, receiver_coords=pillar_1A_base, L1=True)
G17_rover_var = Variance(sat_coords=G17, receiver_coords=pillar_3A_rover,  L1=True)
G15_base_var = Variance(sat_coords=G15, receiver_coords=pillar_1A_base,  L1=True)
G15_rover_var = Variance(sat_coords=G15, receiver_coords=pillar_3A_rover, L1=True)
G13_base_var = Variance(sat_coords=G13, receiver_coords=pillar_1A_base,  L1=True)
G13_rover_var = Variance(sat_coords=G13, receiver_coords=pillar_3A_rover, L1=True)
G12_base_var = Variance(sat_coords=G12, receiver_coords=pillar_1A_base,  L1=True)
G12_rover_var = Variance(sat_coords=G12, receiver_coords=pillar_3A_rover, L1=True)
G10_base_var = Variance(sat_coords=G10, receiver_coords=pillar_1A_base, L1=True)
G10_rover_var = Variance(sat_coords=G10, receiver_coords=pillar_3A_rover, L1=True)

variance_vector = np.array([ [G24_base_var .elevation_variance_calculator()],
                             [G24_rover_var.elevation_variance_calculator()],
                             [G19_base_var .elevation_variance_calculator()],
                             [G19_rover_var.elevation_variance_calculator()],
                             [G18_base_var .elevation_variance_calculator()],
                             [G18_rover_var.elevation_variance_calculator()],
                             [G17_base_var .elevation_variance_calculator()],
                             [G17_rover_var.elevation_variance_calculator()],
                             [G15_base_var .elevation_variance_calculator()],
                             [G15_rover_var.elevation_variance_calculator()],
                             [G13_base_var .elevation_variance_calculator()],
                             [G13_rover_var.elevation_variance_calculator()],
                             [G12_base_var .elevation_variance_calculator()],
                             [G12_rover_var.elevation_variance_calculator()],
                             [G10_base_var .elevation_variance_calculator()],
                             [G10_rover_var.elevation_variance_calculator()]])


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


def Cd_calculator(D, S, Cl):
    return (((D.dot(S)).dot(Cl)).dot(transpose(S))).dot(transpose(D))


def calculate_x_hat(A, W, b):
    return ((linalg.inv((transpose(A).dot(W)).dot(A))).dot(transpose(A).dot(W))).dot(b)


def ATWA(A, W):
    return ((transpose(A)).dot(W)).dot(A)

if __name__ == "__main__":

    # Suppress Scientific mode for simplicity
    np.set_printoptions(suppress=True,
                        formatter={'float_kind':'{:16.3f}'.format},
                        linewidth=130)

    vec2 = variance_vector.reshape(variance_vector.shape[0], 1)
    ax = sns.heatmap((vec2),
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False,
                fmt='g')
    plt.title("Vector of variances")
    plt.savefig("Matrices/Vector of variances")
    plt.show()



    """
    l vector - The vector of raw observations
    """

    vec2 = l.reshape(l.shape[0], 1)
    ax = sns.heatmap((vec2),
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False,
                fmt='g')
    plt.title("Vector of observations (l) Matrix")
    plt.savefig("Matrices/Vector of observations (l) Matrix")
    plt.show()

    """
    For purposes of single differencing and double differencing.
    The D and S are created.
    D DIM: 7 x 8
    S DIM: 8 X 16
    """

    ax = sns.heatmap(S,
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False)
    plt.title("Single Differencing (S) Matrix")
    plt.savefig("Matrices/Single Differencing (S) Matrix")
    plt.show()

    sns.heatmap(D,
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False)
    plt.title("Double Differencing (D) Matrix")
    plt.savefig("Matrices/Double Differencing (D) Matrix")
    plt.show()


    """
    SINGLE DIFFERENCING
    As receiver 1A and 3A are relatively close together (9.4m)
    The ionspheric and tropospheric terms cancel out
    Therefore we can calculate receiver-receiver single difference (RRSD)
    sl = vector of receiver-receiver single differences 
    DIM: 8 x 1
    """

    sl = S.dot(l)
    vec2 = sl.reshape(sl.shape[0], 1)
    ax = sns.heatmap((vec2),
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False)
    plt.title("Vector of single differences (sl) Matrix")
    plt.savefig("Matrices/Vector of single differences (sl) Matrix")
    plt.show()

    """
    DOUBLE DIFFERENCING 
    
    We choose the satellite with the highest elevation as the reference satellite.
    This will be G24 - This satellite will have the least noise and multipath 
    
    Dsl = vector of double differences
    DIM: 7 x 1
    """

    Dsl = D.dot(sl)
    vec2 = Dsl.reshape(Dsl.shape[0], 1)
    ax = sns.heatmap((vec2),
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False)
    plt.title("Vector of Double differences (Dsl) Matrix")
    plt.savefig("Matrices/Vector of Double differences (Dsl) Matrix")
    plt.show()


    """
    Covariance matrix of the observation vector.
    This is an identity matrix scaled by the variance.
    It is assumed that the variance of each raw phase observation has the same l1 standard deviation is 0.003. 
    Therefore variance is this value squared.
    DIM: 16 x 16 
    """

    cl = np.eye(16, 16) * variance_vector
    ax = sns.heatmap((cl),
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False)
    plt.title("Covariance matrix of observations (cl) Matrix")
    plt.savefig("Matrices/Covariance matrix of observations (cl) Matrix")
    plt.show()



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

    Cd = (Cd_calculator(D, S, cl))
    ax = sns.heatmap((Cd),
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False)
    plt.title("covariance matrix of the double differences (Cd) Matrix")
    plt.savefig("Matrices/covariance matrix of the double differences (Cd) Matrix")
    plt.show()




    """
    Calculate the weight matrix of the double differences. 
    
    Wd = Cd^-1
    DIM: 7 x 7
    """

    Wd = linalg.inv(Cd)

    Wd = np.around(Wd, decimals=4)
    heat_map = sns.heatmap(Wd,
                           annot=True,
                           xticklabels=False,
                           yticklabels=False,
                           cmap="Blues",
                           cbar=False, )
    plt.title("Weight (Wd) Matrix")
    plt.savefig("Matrices/Weight (Wd) Matrix")
    plt.show()

    """
    3 methods for calculating the partial differentials of the design matrix - (X Y and Z)
    1 method for calculating the b vector i.e The observed - computed measurements. 
    
    ATTRIBUTES: 
    L1: - If L1 is True, then the wavelength for this 
    brrs: Base Receiver to Reference Satellite geometric range (m)
    rrrs: Roving Receiver to Reference Satellite geometric range (m)
    brcs: Base Receiver to Corresponding Satellite  geometric range (m)
    rrcs: Roving Receiver to Corresponding Satellite geometric range (m)
    N: Ambiguity Term - Always an integer
    e: Noise term
    ref_station: Cartesian [X, Y, Z] coordinates <-- Must be in this order
    Corresponding_sat: Cartesian [X, Y, Z] coordinates <-- Must be in this order
    sat_ref: Cartesian [X, Y, Z] coordinates <-- Must be in this order
    Dsl: Vector of the observed double differences 
    """
    G24G19 = DD(L1=True,
                brrs=G24_base_obs[0],
                rrrs=G24_rover_obs[0],
                brcs=G19_base_obs[0],
                rrcs=G19_rover_obs[0],
                N=G24toG19_after,
                e=G24toG19_noise,
                ref_station=pillar_1A_base,
                corresponding_sat=G19,
                sat_ref=G24,
                dsl=Dsl[0])

    G24G18 = DD(L1=True,
                brrs=G24_base_obs[0],
                rrrs=G24_rover_obs[0],
                brcs=G18_base_obs[0],
                rrcs=G18_rover_obs[0],
                N=G24toG18_after,
                e=G24toG18_noise,
                ref_station=pillar_1A_base,
                corresponding_sat=G18,
                sat_ref=G24,
                dsl=Dsl[1])

    G24G17 = DD(L1=True,
                brrs=G24_base_obs[0],
                rrrs=G24_rover_obs[0],
                brcs=G17_base_obs[0],
                rrcs=G17_rover_obs[0],
                N=G24toG17_after,
                e=G24toG17_noise,
                ref_station=pillar_1A_base,
                corresponding_sat=G17,
                sat_ref=G24, dsl=Dsl[2])

    G24G15 = DD(L1=True,
                brrs=G24_base_obs[0],
                rrrs=G24_rover_obs[0],
                brcs=G15_base_obs[0],
                rrcs=G15_rover_obs[0],
                N=G24toG15_after,
                e=G24toG15_noise,
                ref_station=pillar_1A_base,
                corresponding_sat=G15,
                sat_ref=G24, dsl=Dsl[3])

    G24G13 = DD(L1=True,
                brrs=G24_base_obs[0],
                rrrs=G24_rover_obs[0],
                brcs=G13_base_obs[0],
                rrcs=G13_rover_obs[0],
                N=G24toG13_after,
                e=G24toG13_noise,
                ref_station=pillar_1A_base,
                corresponding_sat=G13,
                sat_ref=G24,
                dsl=Dsl[4])

    G24G12 = DD(L1=True,
                brrs=G24_base_obs[0],
                rrrs=G24_rover_obs[0],
                brcs=G12_base_obs[0],
                rrcs=G12_rover_obs[0],
                N=G24toG12_after,
                e=G24toG12_noise,
                ref_station=pillar_1A_base,
                corresponding_sat=G12,
                sat_ref=G24,
                dsl=Dsl[5])

    G24G10 = DD(L1=True,
                brrs=G24_base_obs[0],
                rrrs=G24_rover_obs[0],
                brcs=G10_base_obs[0],
                rrcs=G10_rover_obs[0],
                N=G24toG10_after,
                e=G24toG10_noise,
                ref_station=pillar_1A_base,
                corresponding_sat=G10,
                sat_ref=G24,
                dsl=Dsl[6])


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


    vec2 = b.reshape(Dsl.shape[0], 1)
    ax = sns.heatmap((b),
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False)
    plt.title("Observed - Computed (b) Vector")
    plt.savefig("Matrices/Observed - Computed (b) Vector")
    plt.show()



    A = np.array([[G24G19.x_diff(), G24G19.y_diff(), G24G19.z_diff()],
                  [G24G18.x_diff(), G24G18.y_diff(), G24G18.z_diff()],
                  [G24G17.x_diff(), G24G17.y_diff(), G24G17.z_diff()],
                  [G24G15.x_diff(), G24G15.y_diff(), G24G15.z_diff()],
                  [G24G13.x_diff(), G24G13.y_diff(), G24G13.z_diff()],
                  [G24G12.x_diff(), G24G12.y_diff(), G24G12.z_diff()],
                  [G24G10.x_diff(), G24G10.y_diff(), G24G10.z_diff()]])

    A = np.around(A, decimals=4)
    sns.heatmap(A,
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False)
    plt.title("Design (A) Matrix")
    plt.savefig("Matrices/Design (A) Matrix")
    plt.show()

    """
    Output the ATWA matrix 
    """
    atwa = ATWA(A, Wd)
    sns.heatmap(atwa,
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False)
    plt.title("ATWA Matrix")
    plt.savefig("Matrices/ATWA Matrix")
    plt.show()

    """
    Output the (ATWA)^-1 matrix 
    """
    inverse_atwa = linalg.inv(atwa)
    sns.heatmap(inverse_atwa,
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False)
    plt.title("(ATWA)^-1 Matrix")
    plt.savefig("Matrices/(ATWA)^-1 Matrix")
    plt.show()


    # Calculate X_hat
    X_hat = calculate_x_hat(A, Wd, b)

    L1_wl = 0.19

    X_hat = [(L1_wl * X_hat[0]), (L1_wl * X_hat[1]), (L1_wl * X_hat[2])]

    updated_pillar_3A = [(pillar_3A_rover[0] + X_hat[0]), (pillar_3A_rover[1] + X_hat[1]), (pillar_3A_rover[2] + X_hat[2])]



    """
    Table to display final results
    """
    table = np.array([["X", float(pillar_3A_rover[0]), float(X_hat[0]), float(updated_pillar_3A[0])],
                     ["Y", float(pillar_3A_rover[1]), float(X_hat[1]), float(updated_pillar_3A[1])],
                      ["Z", float(pillar_3A_rover[2]), float(X_hat[2]), float(updated_pillar_3A[2])]])

    fig, ax = plt.subplots()
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')
    columns = [" ", "Nominal Coordinates", "Updates", "Updated Coordinates"]
    df = pd.DataFrame(table, columns=columns)
    ax.table(cellText=df.values, colLabels=df.columns, loc='center')
    fig.tight_layout()
    plt.savefig("Matrices/Nominal and Updated Coords")
    plt.show()


    """
    As a performance test, the computed distances between pillar 1A and pillar 3A are compared.
    """
    print(distance(after_ambiguity_resolution, pillar_1A_base))

    table = np.array([["Nominal: ", distance(pillar_1A_base, pillar_3A_rover)],
                      ["Updated: ", distance(pillar_1A_base, updated_pillar_3A)]])
    fig, ax = plt.subplots()
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')
    columns = [" ", "Distance between Pillars (m)"]
    df = pd.DataFrame(table, columns=columns)
    ax.table(cellText=df.values, colLabels=df.columns, loc='center')
    fig.tight_layout()
    plt.savefig("Matrices/Nominal Vs Updated Distances")
    plt.show()


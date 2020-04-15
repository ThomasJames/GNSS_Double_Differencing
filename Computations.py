from math import sqrt, cos, sin, degrees, acos
import numpy as np
from numpy import transpose, linalg
from Matrix_Computation_Classes import DD, HeatMap, MatrixOperations

print(np.eye(16, 16))        

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
before_pa =       np.array([[4929605.096], [-29123.627], [4033604.055]])
after_pa =        np.array([[4929604.918], [-29123.616], [4033603.990]])
distance_between_receivers = 94.4  # Measured in meters / Approximate

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


def elevation_variance_calculator(sat_coords, receiver_coords, range_obs):

    # Extract X, Y, Z
    x_s, y_s, z_s = sat_coords[0], sat_coords[1], sat_coords[2]
    # x_r, y_r, z_r = receiver_coords[0], receiver_coords[1], receiver_coords[2]

    # Calculate distance from EC
    ec_s = sqrt((sqrt(x_s**2 + y_s**2))**2 + z_s**2)
    # ec_r = sqrt((sqrt(x_s ** 2 + y_s ** 2))**2 + z_r**2)

    # Using known approximation of earth radius.
    ec_r = 6378000

    # Calculate the earth surface angle
    angle = degrees(acos((ec_r**2 + range_obs[0]**2 - ec_s**2) / (2 * ec_r * range_obs[0])))

    # L1 standard deviation
    l1_SD = 0.003

    # Calculate variance 
    variance = (l1_SD ** 2) / sin(angle)

    return variance

G24_variance_base = elevation_variance_calculator(G24, pillar_1A_base, G24_base_obs)
G24_variance_rover = elevation_variance_calculator(G24, pillar_3A_rover, G24_rover_obs)
G19_variance_base = elevation_variance_calculator(G19, pillar_1A_base, G19_base_obs)
G19_variance_rover = elevation_variance_calculator(G19, pillar_3A_rover, G19_rover_obs)
G18_variance_base = elevation_variance_calculator(G18, pillar_1A_base, G18_base_obs)
G18_variance_rover = elevation_variance_calculator(G18, pillar_3A_rover, G18_rover_obs)
G17_variance_base = elevation_variance_calculator(G17, pillar_1A_base, G17_base_obs)
G17_variance_rover = elevation_variance_calculator(G17, pillar_3A_rover, G17_rover_obs)
G15_variance_base = elevation_variance_calculator(G15, pillar_1A_base, G15_base_obs)
G15_variance_rover = elevation_variance_calculator(G15, pillar_3A_rover, G15_rover_obs)
G13_variance_base = elevation_variance_calculator(G13, pillar_1A_base, G13_base_obs)
G13_variance_rover = elevation_variance_calculator(G13, pillar_3A_rover, G13_rover_obs)
G12_variance_base = elevation_variance_calculator(G12, pillar_1A_base, G12_base_obs)
G12_variance_rover = elevation_variance_calculator(G12, pillar_3A_rover, G12_rover_obs)
G10_variance_base = elevation_variance_calculator(G10, pillar_1A_base, G10_base_obs)
G10_variance_rover = elevation_variance_calculator(G10, pillar_3A_rover, G10_rover_obs)

variance_vector = np.array([ [G24_variance_base ],
                             [G24_variance_rover],
                             [G19_variance_base ],
                             [G19_variance_rover],
                             [G18_variance_base ],
                             [G18_variance_rover],
                             [G17_variance_base ],
                             [G17_variance_rover],
                             [G15_variance_base ],
                             [G15_variance_rover],
                             [G13_variance_base ],
                             [G13_variance_rover],
                             [G12_variance_base ],
                             [G12_variance_rover],
                             [G10_variance_base ],
                             [G10_variance_rover]])


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

    """
    For purposes of single differencing and double differencing.
    The D and S are created.
    D DIM: 7 x 8
    S DIM: 8 X 16
    """

    D_out = HeatMap(matrix=D, title="D_matrix")
    # D_out.output_png()

    S_out = HeatMap(matrix=S, title="S_matrix")
    # S_out.output_png()

    """
    SINGLE DIFFERENCING
    As receiver 1A and 3A are relatively close together (9.4m)
    The ionspheric and tropospheric terms cancel out
    Therefore we can calculate receiver-receiver single difference (RRSD)
    sl = vector of receiver-receiver single differences 
    DIM: 8 x 1
    """

    sl = S.dot(l)

    """
    DOUBLE DIFFERENCING 
    
    We choose the satellite with the highest elevation as the reference satellite.
    This will be G24 - This satellite will have the least noise and multipath 
    
    Dsl = vector of double differences
    DIM: 7 x 1
    """

    Dsl = D.dot(sl)

    """
    Covariance matrix of the observation vector.
    This is an identity matrix scaled by the variance.
    It is assumed that the variance of each raw phase observation has the same l1 standard deviation is 0.003. 
    Therefore variance is this value squared.
    DIM: 16 x 16 
    """

    cl = np.eye(16, 16) * variance_vector
    cl_out = HeatMap(matrix=cl, title="cl_Matrix")
    # cl_out.output_png()

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
    cd_out = HeatMap(matrix=Cd, title="Cd_Matrix")
    # cd_out.output_png()

    """
    Calculate the weight matrix of the double differences. 
    
    Wd = Cd^-1
    DIM: 7 x 7
    """

    Wd = linalg.inv(Cd)
    Wd_out = HeatMap(matrix=Wd, title="Wd_Matrix")
    # Wd_out.output_png()

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

    A = np.array([[G24G19.x_diff(), G24G19.y_diff(), G24G19.z_diff()],
                  [G24G18.x_diff(), G24G18.y_diff(), G24G18.z_diff()],
                  [G24G17.x_diff(), G24G17.y_diff(), G24G17.z_diff()],
                  [G24G15.x_diff(), G24G15.y_diff(), G24G15.z_diff()],
                  [G24G13.x_diff(), G24G13.y_diff(), G24G13.z_diff()],
                  [G24G12.x_diff(), G24G12.y_diff(), G24G12.z_diff()],
                  [G24G10.x_diff(), G24G10.y_diff(), G24G10.z_diff()]])
    # Calculate X_hat

    X_hat = calculate_x_hat(A, Wd, b)

    print("X: ", X_hat[0], "Y:",  X_hat[1], "Z:", X_hat[2])

    print(before_ambiguity_resolution[0] + X_hat[0])
    print(before_ambiguity_resolution[1] + X_hat[1])
    print(before_ambiguity_resolution[2] + X_hat[2])








    


    

    

    










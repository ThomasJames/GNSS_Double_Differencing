from Data import *

"""
GOAL: Calculate the coordinates of the reference antenna (ARP) of the roving receiver 

Try to get your answer close to the figures for 3A. The nominal coordinates given mean you do not need to iterate the 
least squares solution, you should converge on the answer with on round of matrix inversion
"""

def x_differential(reference_station, satelite_corresponding, satelite_reference, wavelength):

    # Condense Coordinates
    X_3A = reference_station[0]
    Y_3A = reference_station[1]
    Z_3A = reference_station[2]

    X_s = satelite_corresponding[0]
    Y_s = satelite_corresponding[1]
    Z_s = satelite_corresponding[2]

    X_s_ref = satelite_reference[0]
    Y_s_ref = satelite_reference[1]
    Z_s_ref = satelite_reference[2]

    result = 1 / wavelength * \
             (
                     (X_3A - X_s) /
                     (sqrt((X_s - X_3A) ** 2 + (Y_s - Y_3A) ** 2 + (Z_s - Z_3A) ** 2))
                     -
                     (X_3A - X_s_ref) /
                     (sqrt(X_s_ref - X_3A) ** 2 + (Y_s_ref - Y_3A) ** 2 + (Z_s_ref - Z_3A) ** 2)
             )
    return float(result)


def y_differential(reference_station,
                   satelite_corresponding,
                   satelite_reference,
                   wavelength):

    # Condense Variables

    X_3A = reference_station[0]
    Y_3A = reference_station[1]
    Z_3A = reference_station[2]

    X_s = satelite_corresponding[0]
    Y_s = satelite_corresponding[1]
    Z_s = satelite_corresponding[2]

    X_s_ref = satelite_reference[0]
    Y_s_ref = satelite_reference[1]
    Z_s_ref = satelite_reference[2]

    # Calculate
    result = 1 / wavelength * \
             (
                     (Y_3A - Y_s) /
                     (sqrt((X_s - X_3A) ** 2 + (Y_s - Y_3A) ** 2 + (Z_s - Z_3A) ** 2))
                     -
                     (Y_3A - Y_s_ref) /
                     (sqrt(X_s_ref - X_3A) ** 2 + (Y_s_ref - Y_3A) ** 2 + (Z_s_ref - Z_3A) ** 2)
             )
    return float(result)


def z_differential(reference_station,
                   satelite_corresponding,
                   satelite_reference,
                   wavelength):

    # Condense variables
    X_3A = reference_station[0]
    Y_3A = reference_station[1]
    Z_3A = reference_station[2]

    X_s = satelite_corresponding[0]
    Y_s = satelite_corresponding[1]
    Z_s = satelite_corresponding[2]

    X_s_ref = satelite_reference[0]
    Y_s_ref = satelite_reference[1]
    Z_s_ref = satelite_reference[2]

    result = 1 / wavelength * \
             (
                     (Z_3A - Z_s) /
                     (sqrt((X_s - X_3A) ** 2 + (Y_s - Y_3A) ** 2 + (X_s - X_3A) ** 2 + (Z_s - Z_3A) ** 2))
                     -
                     (Z_3A - Z_s_ref) /
                     (sqrt(X_s_ref - X_3A) ** 2 + (Y_s_ref - Y_3A) ** 2 + (Z_s_ref - Z_3A) ** 2)
             )
    return float(result)


"""
b_vector function
"""
def b_vector(
        base_range_ref,
        base_range_corresponding,
        rover_range_corresponding,
        rover_range_ref,
        pa,
        wl,
        obs):

    # Condense variables
    brf = base_range_ref[0]
    brc = base_range_corresponding[0]
    rrr = rover_range_ref[0]
    rrc = rover_range_corresponding[0]

    result = obs - 1 / wl * (brf - rrr - brc + rrc) - pa

    return result


def Cd_calculator(D, S, Cl):

    result = (((D.dot(S)).dot(Cl)).dot(transpose(S))).dot(transpose(D))
    return result



def calculate_x_hat(A, W, b):
    result = ((linalg.inv((transpose(A).dot(W)).dot(A))).dot(transpose(A).dot(W))).dot(b)
    return result






if __name__ == "__main__":
    """    
    We will carry out the double differencing with respect to one satellite.
    We use the one with the highest elevation, as this satellite will pick up the least noise:
    1. Short atmospheric path length. 
    2. Maximum reduction of multipath.
    
    Typical observation rates might be every second / 5s / 10s / 30s.
    2 stations 
    8 satellites 
    60 seconds per minute 
    If the A to B was tracked for 10 minutes: 
    2 x 8 x 60 x 10 = 9600
    9600 measurements
    """

    print("vector of observations: ", l)
    print("  ")

    # Calculate the vector of single differences
    sl = S.dot(l)

    # Calculate the vector of double differences
    Dsl = D.dot(sl)


    # Constructing the Design Matrix
    A = np.array([

    [x_differential(pillar_1A_base, G10, G24, wavelength),
     y_differential(pillar_1A_base, G10, G24, wavelength),
     z_differential(pillar_1A_base, G10, G24, wavelength), 1, 0, 0, 0, 0, 0, 0],

    [x_differential(pillar_1A_base, G12, G24, wavelength),
     y_differential(pillar_1A_base, G12, G24, wavelength),
     z_differential(pillar_1A_base, G12, G24, wavelength), 0, 1, 0, 0, 0, 0, 0],

    [x_differential(pillar_1A_base, G13, G24, wavelength),
     y_differential(pillar_1A_base, G13, G24, wavelength),
     z_differential(pillar_1A_base, G13, G24, wavelength), 0, 0, 1, 0, 0, 0, 0],

    [x_differential(pillar_1A_base, G15, G24, wavelength),
     y_differential(pillar_1A_base, G15, G24, wavelength),
     z_differential(pillar_1A_base, G15, G24, wavelength), 0, 0, 0, 1, 0, 0, 0],

    [x_differential(pillar_1A_base, G17, G24, wavelength),
     y_differential(pillar_1A_base, G17, G24, wavelength),
     z_differential(pillar_1A_base, G17, G24, wavelength), 0, 0, 0, 0, 1, 0, 0],

    [x_differential(pillar_1A_base, G18, G24, wavelength),
     y_differential(pillar_1A_base, G18, G24, wavelength),
     z_differential(pillar_1A_base, G18, G24, wavelength), 0, 0, 0, 0, 0, 1, 0],

    [x_differential(pillar_1A_base, G19, G24, wavelength),
     y_differential(pillar_1A_base, G19, G24, wavelength),
     z_differential(pillar_1A_base, G19, G24, wavelength), 0, 0, 0, 0, 0, 0, 1],
    ])

    print("The A matrix: ", A)
    print("  ")

    cl = l1c_variance * np.eye(16 , 16)
    print("The covariance matrix is: ", cl) # Exists and is known
    print("  ")

    print("D:", D.shape)
    print("  ")
    print("S:", S.shape)
    print("  ")
    print("cl: ", cl.shape)
    print("  ")

    Cd = (Cd_calculator(D, S, cl))
    print(Cd)
    print(Cd.shape)
    print(" ")

    Wd = linalg.inv(Cd)
    print("Weight Matrix d", Wd)

    calculate_x_hat(A, Wd, b)












    











    











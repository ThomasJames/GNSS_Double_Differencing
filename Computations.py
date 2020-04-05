from Data import *

a = 1

"""
GOAL: Calculate the coordinates of the reference antenna (ARP) of the roving receiver 

Try to get your answer close to the figures for 3A. The nominal coordinates given mean you do not need to iterate the 
least squares solution, you should converge on the answer with on round of matrix inversion
"""

def x_differential(reference_station, satelite_corresponding, satelite_reference, wavelength):

    # Condense Variables
    X_3A, Y_3A, Z_3A = reference_station[0], reference_station[1], reference_station[2]
    X_s, Y_s, Z_s = satelite_corresponding[0], satelite_corresponding[1], satelite_corresponding[2]
    X_s_ref, Y_s_ref, Z_s_ref = satelite_reference[0], satelite_reference[1], satelite_reference[2]

    result = 1 / wavelength * \
             (
                     (X_3A - X_s) /
                     (sqrt((X_s - X_3A) ** 2 + (Y_s - Y_3A) ** 2 + (Z_s - Z_3A) ** 2))
                     -
                     (X_3A - X_s_ref) /
                     (sqrt(X_s_ref - X_3A) ** 2 + (Y_s_ref - Y_3A) ** 2 + (Z_s_ref - Z_3A) ** 2)
             )
    return float(result)


def y_differential(reference_station, satelite_corresponding, satelite_reference, wavelength):

    # Condense Variables
    # Condense Variables
    X_3A, Y_3A, Z_3A = reference_station[0], reference_station[1], reference_station[2]
    X_s, Y_s, Z_s = satelite_corresponding[0], satelite_corresponding[1], satelite_corresponding[2]
    X_s_ref, Y_s_ref, Z_s_ref = satelite_reference[0], satelite_reference[1], satelite_reference[2]

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


def z_differential(reference_station, satelite_corresponding, satelite_reference, wavelength):

    # Condense Variables
    # Condense Variables
    X_3A, Y_3A, Z_3A = reference_station[0], reference_station[1], reference_station[2]
    X_s, Y_s, Z_s = satelite_corresponding[0], satelite_corresponding[1], satelite_corresponding[2]
    X_s_ref, Y_s_ref, Z_s_ref = satelite_reference[0], satelite_reference[1], satelite_reference[2]

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
pa: Phase ambiguity 
obs: G24 to G10 after ambiguity resolution 
Wl: defined 
"""
def b_vector( base_range_ref,
              base_range_corresponding,
              rover_range_ref,
              rover_range_corresponding,
              N,
              wl,
              obs):

    # Condense variables
    brf = base_range_ref[0]
    brc = base_range_corresponding[0]
    rrr = rover_range_ref[0]
    rrc = rover_range_corresponding[0]

    result = obs - (1 / wl * (brf - rrr - brc + rrc) - N)
    return result

def Cd_calculator(D, S, Cl):
    result = (((D.dot(S)).dot(Cl)).dot(transpose(S))).dot(transpose(D))
    return result

def calculate_x_hat(A, W, b):
    result = ((linalg.inv((transpose(A).dot(W)).dot(A))).dot(transpose(A).dot(W))).dot(b)
    return result

def calculate_measured(wavelength, br, rr, bc, rc, pa, noise):
    cbr = br[1]
    crr = rr[1]
    cbc = bc[1]
    crc = rc[1]

    result = (1 / wavelength * (cbr - crr - cbc + crc)) + pa + noise
    return result
    

if __name__ == "__main__":


    print("vector of observations: ", l)
    print("  ")

    # Calculate the vector of single differences
    sl = S.dot(l)

    print("sl:", sl)

    # Calculate the vector of double differences
    Dsl = D.dot(sl)
    print("Dsl: ", Dsl)
                                 
    # Constructing the Design Matrix
    A = np.array([

    [x_differential(pillar_1A_base, G19, G24, wl),
     y_differential(pillar_1A_base, G19, G24, wl),
     z_differential(pillar_1A_base, G19, G24, wl), 1, 0, 0, 0, 0, 0, 0],

    [x_differential(pillar_1A_base, G18, G24, wl),
     y_differential(pillar_1A_base, G18, G24, wl),
     z_differential(pillar_1A_base, G18, G24, wl), 0, 1, 0, 0, 0, 0, 0],

    [x_differential(pillar_1A_base, G17, G24, wl),
     y_differential(pillar_1A_base, G17, G24, wl),
     z_differential(pillar_1A_base, G17, G24, wl), 0, 0, 1, 0, 0, 0, 0],

    [x_differential(pillar_1A_base, G15, G24, wl),
     y_differential(pillar_1A_base, G15, G24, wl),
     z_differential(pillar_1A_base, G15, G24, wl), 0, 0, 0, 1, 0, 0, 0],

    [x_differential(pillar_1A_base, G13, G24, wl),
     y_differential(pillar_1A_base, G13, G24, wl),
     z_differential(pillar_1A_base, G13, G24, wl), 0, 0, 0, 0, 1, 0, 0],

    [x_differential(pillar_1A_base, G12, G24, wl),
     y_differential(pillar_1A_base, G12, G24, wl),
     z_differential(pillar_1A_base, G12, G24, wl), 0, 0, 0, 0, 0, 1, 0],

    [x_differential(pillar_1A_base, G10, G24, wl),
     y_differential(pillar_1A_base, G10, G24, wl),
     z_differential(pillar_1A_base, G10, G24, wl), 0, 0, 0, 0, 0, 0, 1],
    ])

    # print("A :", A)
    print("A shape:", A.shape)
    print("AT shape: ", transpose(A).shape)
    # print("A inverse shape", (np.linalg.inv(A)).shape)

    # Calculate the observations
    G24toG10_measured = calculate_measured(wl, G24_base_obs, G24_rover_obs, G19_base_obs,
                                           G19_rover_obs, G24toG19_after, G24toG19_noise)
    G24toG12_measured = calculate_measured(wl, G24_base_obs, G24_rover_obs, G18_base_obs,
                                           G18_rover_obs, G24toG18_after, G24toG18_noise)
    G24toG13_measured = calculate_measured(wl, G24_base_obs, G24_rover_obs, G17_base_obs,
                                           G17_rover_obs, G24toG17_after, G24toG17_noise)
    G24toG15_measured = calculate_measured(wl, G24_base_obs, G24_rover_obs, G15_base_obs,
                                           G15_rover_obs, G24toG15_after, G24toG15_noise)
    G24toG17_measured = calculate_measured(wl, G24_base_obs, G24_rover_obs, G13_base_obs,
                                           G13_rover_obs, G24toG13_after, G24toG13_noise)
    G24toG18_measured = calculate_measured(wl, G24_base_obs, G24_rover_obs, G12_base_obs,
                                           G12_rover_obs, G24toG12_after, G24toG12_noise)
    G24toG19_measured = calculate_measured(wl, G24_base_obs, G24_rover_obs, G10_base_obs,
                                           G10_rover_obs, G24toG10_after, G24toG10_noise)




    b = np.array([
        [b_vector(G24_base_obs, G19_base_obs, G24_rover_obs, G19_rover_obs,  G24toG19_before, wl, G24toG19_measured)],
        [b_vector(G24_base_obs, G18_base_obs, G24_rover_obs, G18_rover_obs,  G24toG18_before, wl, G24toG18_measured)],
        [b_vector(G24_base_obs, G17_base_obs, G24_rover_obs, G17_rover_obs,  G24toG17_before, wl, G24toG17_measured)],
        [b_vector(G24_base_obs, G15_base_obs, G24_rover_obs, G15_rover_obs,  G24toG15_before, wl, G24toG15_measured)],
        [b_vector(G24_base_obs, G13_base_obs, G24_rover_obs, G13_rover_obs,  G24toG13_before, wl, G24toG13_measured)],
        [b_vector(G24_base_obs, G12_base_obs, G24_rover_obs, G12_rover_obs,  G24toG12_before, wl, G24toG12_measured)],
        [b_vector(G24_base_obs, G10_base_obs, G24_rover_obs, G10_rover_obs,  G24toG10_before, wl, G24toG10_measured)],
    ])


    print("b :", b)
    print("b dimensions: ", b.shape)

    # print("The A matrix: ", A)
    print("A dimensions", A.shape)
    print("  ")


    # print("The covariance matrix is: ", cl) # Exists and is known
    print("  ")
    print("D dimensions: ", D.shape)
    print("  ")
    print("S dimensions: ", S.shape)
    print("  ")

    # Calculate the covariance matrix of the observation vector l
    cl = variances * np.eye(8 , 8)
    print("cl", cl)
    print("cl dimensions: ", cl.shape)
    print("  ")
    plt.imshow(cl)
    plt.title('Cov l')
    plt.show()

    Cd = (Cd_calculator(D, S, cl))
    # print(Cd)
    print(Cd.shape)
    print(" ")
    plt.imshow(Cd)
    plt.title('Cov d')
    plt.show()

    Wd = linalg.inv(Cd)
    print("Weight Matrix d", Wd)
    print("Weight shape: ", Wd.shape)
    print("Weight T shape: ", (transpose(Wd)).shape)


    sns.heatmap(Wd, annot=True, fmt="d")
    plt.show()

    x_hat = calculate_x_hat(A, Wd, b)
    # print("x_hat: ", x_hat)

    X = pillar_3A_rover[0] + x_hat[0]
    Y = pillar_3A_rover[1] + x_hat[1]
    Z = pillar_3A_rover[2] + x_hat[2]

    print(X, Y, Z)






































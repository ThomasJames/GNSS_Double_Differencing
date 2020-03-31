import math
import numpy as np
from Data import *

"""
Try to get your answer close to the figures for 3A. The nominal coordinates given mean you do not need to iterate the 
least squares solution, you should converge on the answer with on round of matrix inversion
"""

# 16 x 8:  Differencing matrix
s = np.array([[1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
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


def x_differential(reference_station, satelite_corresponding, satelite_reference, wavelenth):
    # Extract Coordinates
    X_3A = reference_station[0]
    Y_3A = reference_station[1]
    Z_3A = reference_station[2]

    X_s = satelite_corresponding[0]
    Y_s = satelite_corresponding[1]
    Z_s = satelite_corresponding[2]

    X_s_ref = satelite_reference[0]
    Y_s_ref = satelite_reference[1]
    Z_s_ref = satelite_reference[2]

    result = 1 / wavelenth * \
             (
                     (X_3A - X_s) /
                     (math.sqrt((X_s - X_3A) ** 2 + (Y_s - Y_3A) ** 2 + (X_s - X_3A) ** 2 + (Z_s - Z_3A) ** 2))
                     -
                     (X_3A - X_s_ref) /
                     (math.sqrt(X_s_ref - X_3A) ** 2 + (Y_s_ref - Y_3A) ** 2 + (X_s_ref - X_3A) ** 2 +
                      (Z_s_ref - Z_3A) ** 2)
             )
    return float(result)


def y_differential(reference_station, satelite_corresponding, satelite_reference, wavelenth):
    # Extract Coordinates
    X_3A = reference_station[0]
    Y_3A = reference_station[1]
    Z_3A = reference_station[2]

    X_s = satelite_corresponding[0]
    Y_s = satelite_corresponding[1]
    Z_s = satelite_corresponding[2]

    X_s_ref = satelite_reference[0]
    Y_s_ref = satelite_reference[1]
    Z_s_ref = satelite_reference[2]

    result = 1 / wavelenth * \
             (
                     (Y_3A - Y_s) /
                     (math.sqrt((X_s - X_3A) ** 2 + (Y_s - Y_3A) ** 2 + (X_s - X_3A) ** 2 + (Z_s - Z_3A) ** 2))
                     -
                     (Y_3A - Y_s_ref) /
                     (math.sqrt(X_s_ref - X_3A) ** 2 + (Y_s_ref - Y_3A) ** 2 + (X_s_ref - X_3A) ** 2 +
                      (Z_s_ref - Z_3A) ** 2)
             )
    return float(result)


def z_differential(reference_station, satelite_corresponding, satelite_reference, wavelenth):
    # Extract Coordinates
    X_3A = reference_station[0]
    Y_3A = reference_station[1]
    Z_3A = reference_station[2]

    X_s = satelite_corresponding[0]
    Y_s = satelite_corresponding[1]
    Z_s = satelite_corresponding[2]

    X_s_ref = satelite_reference[0]
    Y_s_ref = satelite_reference[1]
    Z_s_ref = satelite_reference[2]

    result = 1 / wavelenth * \
             (
                     (Z_3A - Z_s) /
                     (math.sqrt((X_s - X_3A) ** 2 + (Y_s - Y_3A) ** 2 + (X_s - X_3A) ** 2 + (Z_s - Z_3A) ** 2))
                     -
                     (Z_3A - Z_s_ref) /
                     (math.sqrt(X_s_ref - X_3A) ** 2 + (Y_s_ref - Y_3A) ** 2 + (X_s_ref - X_3A) ** 2 +
                      (Z_s_ref - Z_3A) ** 2)
             )
    return float(result)


if __name__ == "__main__":
    print("Script start")

    """
    Typical observation rates might be every second / 5s / 10s / 30s
    2 stations 
    8 satellites 
    60 seconds per minute 
    If the A to B was tracked for 10 minutes: 
    2 x 8 x 60 x 10 = 9600
    9600 measurements
    We are assuming tha
    """

    print(l)

    # Calculate the vector of single differences
    sl = s.dot(l)

    # Calculate the vector of double differences
    Dsl = D.dot(sl)
    print(Dsl)

    wavelength = G10[1] / G10[0]

    # Constructing the Design Matrix
    design = np.array([

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





    











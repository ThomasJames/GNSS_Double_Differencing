from math import sqrt
from numpy import transpose, linalg
import matplotlib.pyplot as plt
import seaborn as sns
from math import sqrt, cos, sin, degrees, acos, radians
from typing import List, Tuple

"""
Class DD

ref_station - Earth centric Cartesian Coordinates of the reference station [X, Y, Z]
corresponding_sat - Earth centric cartesian Coordinates of the corresponding station [X, Y, Z]
sat_ref - Satellite reference of cartesian Coordinates of the reference station [X, Y, Z]
wl - Wavelength 
brrs - Base receiver to reference satellite 
rrrs - Reference receiver to reference satellite
brcs - Base receiver to corresponding satellite 
rrcs - Reference receiver to corresponding satellite
N - Ambiguity term (int)
e - Noise term 
c = speed of light in vacuum (299792458.0 ms-1) - Set to default    
f = signal frequency (L1: 1575.42MHz, L2: 1227.6MHz)  either L1 or L2 can be True
λ=𝑐/𝑓  - Wavelength calculated from c and f 

"""

class DD:
    def __init__(self, ref_station =None,
                       corresponding_sat =None,
                       sat_ref= None,
                       N: int = None,
                       e: float = None,
                       brrs: float = None,
                       rrrs: float = None,
                       brcs: float = None,
                       rrcs: float =None,
                       L1 : bool = True,
                       dsl: float = None):

        if L1:
            wl = 0.19029367


        self.X_3A = ref_station[0]
        self.Y_3A = ref_station[1]
        self.Z_3A = ref_station[2]
        self.X_s = corresponding_sat[0]
        self.Y_s = corresponding_sat[1]
        self.Z_s = corresponding_sat[2]
        self.X_s_ref = sat_ref[0]
        self.Y_s_ref = sat_ref[1]
        self.Z_s_ref = sat_ref[2]
        self.wl = wl
        self.N = N
        self.e = e
        self.brrs = brrs
        self.rrrs = rrrs
        self.brcs = brcs
        self.rrcs = rrcs
        self.dsl = dsl


    def x_diff(self) -> float:
        return float(1 / self.wl * \
                    (
                         (self.X_3A - self.X_s) /
                         (sqrt((self.X_s - self.X_3A) ** 2 +
                               (self.Y_s - self.Y_3A) ** 2 +
                               (self.Z_s - self.Z_3A) ** 2))
                         -
                         (self.X_3A - self.X_s_ref) /
                         (sqrt(self.X_s_ref - self.X_3A) ** 2 +
                          (self.Y_s_ref - self.Y_3A) ** 2 +
                          (self.Z_s_ref - self.Z_3A) ** 2)))

    def y_diff(self) -> float:
        return float((1 / self.wl * \
                    (
                         (self.Y_3A - self.Y_s) /
                         (sqrt((self.X_s - self.X_3A) ** 2 +
                               (self.Y_s - self.Y_3A) ** 2 +
                               (self.Z_s - self.Z_3A) ** 2))
                         -
                         (self.Y_3A - self.Y_s_ref) /
                         (sqrt(self.X_s_ref - self.X_3A) ** 2 +
                          (self.Y_s_ref - self.Y_3A) ** 2 +
                          (self.Z_s_ref - self.Z_3A) ** 2))))

    def z_diff(self) -> float:
        return float(1 / self.wl * \
                     (
                         (self.Z_3A - self.Z_s) /
                         (sqrt((self.X_s - self.X_3A) ** 2 +
                               (self.Y_s - self.Y_3A) ** 2 +
                               (self.Z_s - self.Z_3A) ** 2))
                         -
                         (self.Z_3A - self.Z_s_ref) /
                         (sqrt(self.X_s_ref - self.X_3A) ** 2 +
                          (self.Y_s_ref - self.Y_3A) ** 2 +
                          (self.Z_s_ref - self.Z_3A) ** 2)))

    def calc_b_vector(self) -> float:
        # observed - The vector of measured quantities
        o = self.dsl - self.N + self.e

        # Computed
        c = (1 / self.wl) * (self.brrs - self.rrrs - self.brcs + self.rrcs)
        return o - c


class Variance:

    def __init__(self, sat_coords: List[float],
                       receiver_coords: List[float],
                       L1: bool == True):

        if L1:
            l1_std = 0.003

        assert len(sat_coords) == len(receiver_coords)


        self.l1_std = l1_std
        self.sat_coords = sat_coords
        self.receiver_coords = receiver_coords

    def elevation_variance_calculator(self) -> float:

        """"
        This method calculates the satellite angle of elevation in the following stages:
        Calculates the distance of receiver to the satellite (m) using pythagoras theorem.
        Calculates the distance between the earth center and the satellite (m) using pythagoras theorem.
        Calculates the distance between the earth center and the receiver (m) using pythagoras theorem.
        These ranges make up a scalene triangle, where all ranges are known.
        The low of cosines is used to calculate the angle about the receiver in degrees. 
        90 is subtracted from this angle to get the local elevation angle.     
        The method then calculates the variance of the satellite at the calculated angle.
        returns the variance as a float
        """""

        # Extract ECEF distances (m)
        x_s, y_s, z_s = self.sat_coords[0], self.sat_coords[1], self.sat_coords[2]
        x_r, y_r, z_r = self.receiver_coords[0], self.receiver_coords[1], self.receiver_coords[2]

        # Distance from receiver to satellite (m)
        r_s = sqrt((x_s - x_r)**2 + (y_s - y_r)**2 + (z_s - z_r)**2)

        # Distance from earth center to satellite (m)
        ec_s = sqrt((sqrt(x_s ** 2 + y_s ** 2)) ** 2 + z_s ** 2)

        # Distance from earth center to receiver (m)
        ec_r = sqrt((sqrt(x_r ** 2 + y_r ** 2)) ** 2 + z_r ** 2)

        # Angle from the local horizontal to the satellite (m)
        angle = (degrees(acos((ec_r ** 2 + r_s ** 2 - ec_s ** 2) / (2 * ec_r * r_s)))) - 90
        angle = radians(angle)

        # Variance (uncertainty associated with the satellite) (m)
        variance = (self.l1_std ** 2) / (sin(angle))

        return variance




"""
Matrix Operations

D - Double differencing matrix 
S - Single differencing matrix
Cl - Covariance matrix of observations
A - Design matric
W - Weight matrix
b - B (innovation vector?)
"""

class MatrixOperations:
    def __init__(self, D=None, S=None, Cl=None, A=None, W=None, b=None):
        self.D = D
        self.S = S
        self.Cl = Cl
        self.A = A
        self.W = W
        self.b = b

    def Cd_calculator(self):
        try:
            return (((self.D.dot(self.S)).dot(self.Cl)).dot(transpose(self.S))).dot(transpose(self.D))
        except IOError:
            print("Cd_calculator failed")

    def calculate_x_hat(self):
        try:
            return \
                ((linalg.inv((transpose(self.A).dot(self.W)).dot(self.A))).dot(transpose(self.A).dot(self.W))).dot(self.b)
        except IOError:
            print("Calculate_x_hat failed")

    def ATWA(self):
        try:
            return ((transpose(self.A)).dot(self.W)).dot(self.A)
        except IOError:
            print("ATWA failed")

"""
Heatmap Matrix Plotter)

matrix - Numpy file.
title - Title of the matrix (string)
file_type - Set to png as default.
colour - Set to "Pastel1" as default.
"""

class HeatMap:

    def __init__(self, matrix, filename):
        self.matrix = matrix
        self.filename = filename

    def png_out(self):
        heat_map = sns.heatmap(self.matrix,
                               annot=True,
                               xticklabels=False,
                               yticklabels=False,
                               cmap="Blues",
                               cbar=False,)
        plt.imsave('Matrix.png', self.matrix)





from math import sqrt
from numpy import transpose, linalg
import matplotlib.pyplot as plt
import seaborn as sns


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
    def __init__(self, ref_station=None, corresponding_sat=None, sat_ref=None,
                       N=None, e=None,
                       brrs=None, rrrs=None, brcs=None, rrcs=None, L1=False):

        if L1:
            f = 0.19

        self.X_3A = ref_station[0]
        self.Y_3A = ref_station[1]
        self.Z_3A = ref_station[2]
        self.X_s = corresponding_sat[0]
        self.Y_s = corresponding_sat[1]
        self.Z_s = corresponding_sat[2]
        self.X_s_ref = sat_ref[0]
        self.Y_s_ref = sat_ref[1]
        self.Z_s_ref = sat_ref[2]
        self.wl = f
        self.N = N
        self.e = e
        self.brrs = brrs
        self.rrrs = rrrs
        self.brcs = brcs
        self.rrcs = rrcs

    def x_diff(self):
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

    def y_diff(self):
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

    def z_diff(self):
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

    def calc_b_vector(self, dsl):
        # observed - The vector of measured quantities
        o = dsl

        # Computed
        c = (1 / self.wl) * (self.brrs - self.rrrs - self.brcs + self.rrcs) + self.N + self.e
        return o - c


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

    def __init__(self, matrix, title, file_type=".png", colour="Pastel1"):
        self.matrix = matrix
        self.title = title
        self.file_type = file_type
        self.colour = colour

    def output_png(self):
        sns.heatmap(self.matrix,
                    annot=True,
                    cbar=False,
                    xticklabels=False,
                    yticklabels=False,
                    cmap=self.colour)
        plt.title(str(self.title))
        plt.savefig("Matrix_Output/" + str(self.title) + self.file_type)

    def matrix_plotter(self):
        sns.heatmap(self.matrix,
                    annot=True,
                    cbar=False,
                    xticklabels=False,
                    yticklabels=False,
                    cmap=self.colour)
        plt.title(str(self.title))
        plt.show()

from numpy import transpose, linalg
import matplotlib.pyplot as plt
import seaborn as sns
from math import sqrt, sin, degrees, acos, radians
from typing import List

class DD:
    
    """"
    ref_station - Earth centric Cartesian Coordinates of the reference station [X, Y, Z]
    corresponding_sat - Earth centric cartesian Coordinates of the corresponding station [X, Y, Z]
    sat_ref - Satellite reference of cartesian Coordinates of the reference station [X, Y, Z]
    brrs - Base receiver to reference satellite 
    rrrs - Reference receiver to reference satellite
    brcs - Base receiver to corresponding satellite 
    rrcs - Reference receiver to corresponding satellite
    c = speed of light in vacuum (299792458.0 ms-1) - Set to default    
    f = signal frequency (L1: 1575.42MHz, L2: 1227.6MHz)  either L1 or L2 can be True
    λ=𝑐/𝑓  - Wavelength calculated from c and f 
    """""

    def __init__(self, ref_station: List[float] = None,
                       rov_station: List[float] = None,
                       corresponding_sat: List[float] = None,
                       sat_ref: List[float] = None,
                       L1: bool = True,
                       L2: bool = False,
                       observed: float = None):

        # Speed of light m/s
        c: float = 299792458.0

        # Signal frequency of L1 (MHz)
        L1_f: float= 1575.42 * 1000000

        # Signal frequency of L2
        L2_f: float = 1227.6 * 1000000

        # Set to true by default
        if L1:
            wl = c / L1_f

        if L2:
            wl = c / L2_f

        # Error check the arguments
        assert len(ref_station) == len(rov_station) == len(corresponding_sat) == len(sat_ref)
        assert L1 != L2


        # Compute ranges from satellite coordinates
        brrs: float = self.distance(ref_station, sat_ref)
        rrrs: float = self.distance(rov_station, sat_ref)
        brcs: float = self.distance(ref_station, corresponding_sat)
        rrcs: float = self.distance(rov_station, corresponding_sat)

        # Initialise the class variables
        self.x_1A = ref_station[0]
        self.y_1A = ref_station[1]
        self.z_1A = ref_station[2]
        self.x_3A = rov_station[0]
        self.y_3A = rov_station[1]
        self.z_3A = rov_station[2]
        self.x_s = corresponding_sat[0]
        self.y_s = corresponding_sat[1]
        self.z_s = corresponding_sat[2]
        self.x_s_ref = sat_ref[0]
        self.y_s_ref = sat_ref[1]
        self.z_s_ref = sat_ref[2]
        self.wl = wl
        self.brrs = brrs
        self.rrrs = rrrs
        self.brcs = brcs
        self.rrcs = rrcs
        self.observed = observed


    def x_diff(self) -> float:
        return float(1 / self.wl * \
                     (
                             (self.x_3A - self.x_s) /
                             (sqrt(((self.x_s - self.x_3A) ** 2) +
                                   ((self.y_s - self.y_3A) ** 2) +
                                   ((self.z_s - self.z_3A) ** 2)))
                             -
                             (self.x_3A - self.x_s_ref) /
                             (sqrt(((self.x_s_ref - self.x_3A) ** 2) +
                                   ((self.y_s_ref - self.y_3A) ** 2) +
                                   ((self.z_s_ref - self.z_3A) ** 2)))))

    def y_diff(self) -> float:
        return float((1 / self.wl * \
                      (
                              (self.y_3A - self.y_s) /
                              (sqrt(((self.x_s - self.x_3A) ** 2) +
                                    ((self.y_s - self.y_3A) ** 2) +
                                    ((self.z_s - self.z_3A) ** 2)))
                              -
                              (self.y_3A - self.y_s_ref) /
                              (sqrt(((self.x_s_ref - self.x_3A) ** 2) +
                                    ((self.y_s_ref - self.y_3A) ** 2) +
                                    ((self.z_s_ref - self.z_3A) ** 2))))))

    def z_diff(self) -> float:
        return float(1 / self.wl * \
                     (
                             (self.z_3A - self.z_s) /
                             (sqrt(((self.x_s - self.x_3A) ** 2) +
                                   ((self.y_s - self.y_3A) ** 2) +
                                   ((self.z_s - self.z_3A) ** 2)))
                             -
                             (self.z_3A - self.z_s_ref) /
                             (sqrt(((self.x_s_ref - self.x_3A) ** 2) +
                                   ((self.y_s_ref - self.y_3A) ** 2) +
                                   ((self.z_s_ref - self.z_3A) ** 2)))))

    def calc_b_vector(self) -> float:
        # observed - The vector of measured quantities
        o = self.observed

        # Computed
        c = (1 / self.wl) * (self.brrs - self.rrrs - self.brcs + self.rrcs)
        return o - c

    def distance(self, point_1: List[float], point_2: List[float]) -> float:
        """""
        Find the difference between two points given sets of [X, Y, Z] coordinates.
        """""
        return sqrt((point_2[0] - point_1[0]) ** 2 +
                    (point_2[1] - point_1[1]) ** 2 +
                    (point_2[2] - point_1[2]) ** 2)



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

    def elevation_calculator(self) -> float:

        """"
        This method calculates the satellite angle of elevation in the following stages:
        Calculates the distance of receiver to the satellite (m) using pythagoras theorem.
        Calculates the distance between the earth center and the satellite (m) using pythagoras theorem.
        Calculates the distance between the earth center and the receiver (m) using pythagoras theorem.
        These ranges make up a scalene triangle, where all ranges are known.
        The low of cosines is used to calculate the angle about the receiver in degrees. 
        90 is subtracted from this angle to get the local elevation angle.     
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
        angle_radians = radians((degrees(acos((ec_r ** 2 + r_s ** 2 - ec_s ** 2) / (2 * ec_r * r_s)))) - 90)

        return angle_radians
    
    def variance(self):
        """""
        This method then calculates the variance of the satellite at the calculated angle.
        returns the variance as a float
        """""
        # Variance (uncertainty associated with the satellite) (m)
        variance = (self.l1_std ** 2) / (sin(self.elevation_calculator()))
        
        return variance


class MatrixOperations:
    """"
    Matrix Operations
    
    D - Double differencing matrix 
    S - Single differencing matrix
    Cl - Covariance matrix of observations
    A - Design matrix
    W - Weight matrix
    b - B (innovation vector?)
    """""
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


def matrix_heatmap(matrix, name: str) -> None:
    sns.heatmap(matrix,
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False, )
    plt.title(f"{name} Matrix")
    plt.savefig(f"Matrices/{name} Matrix")
    plt.show()

def vector_heatmap(matrix, name: str):
    sns.heatmap(matrix,
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False, )
    plt.title(f"{name} Matrix")
    plt.savefig(f"Vectors/{name} Vector")
    plt.show()

def flipped_vector_heatmap(data, name: str):
    vec2 = data.reshape(data.shape[0], 1)
    ax = sns.heatmap((vec2),
                annot=True,
                xticklabels=False,
                yticklabels=False,
                cmap="Blues",
                cbar=False,
                fmt='g')
    plt.title(f"{name} Vector")
    plt.savefig(f"Vectors/{name} Vector")
    plt.show()




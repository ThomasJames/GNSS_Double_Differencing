# 📡 GNSS Double Differencing  📡

Double Differencing: Main method of high precision commercial positioning. This method requires two high grade GNSS receivers. In this case a reference, and a (static or kinematic) rover receiver.  

## Scenario

Phase and pseudorange observations measured on a calibration baseline in Valencia, Spain. The sensors used
are geodetic quality Leica receivers using choke ring antennas. The receiver on pillar 1A is treated as
the reference receiver, with the following known ECEF coordinates in metres (XYZ)T
for the phase

centre:
4929635.440
 -29041.877
4033567.846

Use the following nominal coordinates for the phase centre of the sensor on pillar 3A:
4929605.400
 -29123.700
4033603.800

Use double differenced phase measurements, from the first epoch of data only (2016 11 15 22 19
5), to compute the precise coordinates of the pillar 3A sensor phase centre

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Diagrams/DD_Diagram_aid.png" width="500">

## Prerequisites 

``` Python 3``` 

The following libraries must be installed:

``` 
pip install numpy 
pip install matplotlib
pip install seaborn 
```

You can find the computations [here](https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Computations.py)

### Clone this repository:

```
$ git clone https://github.com/ThomasJames/GNSS_Double_Differencing
```

## Data

The original [text file](https://github.com/ThomasJames/GNSS_Data_(text).txt) is in this repository, however the data has been organised and processed into [python](https://github.com/ThomasJames/GNSS_Double_Differencing/Data.py) variables. 

## Method 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Diagrams/Method_Flowchart.png" width="600">

Satellite G24 has the highest elevation of apprixmatley 71 degrees. This satellite is used as the reference satellite.

It is important to initially calculate the elevation angles of each satelite. The error of a satellite is inversely proportial to elevation with respect to a local horizon. Low elevation satellites produce less reliable results, and this needs to be taken into account when formulating a weight matrix.


<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Diagrams/Satellite_Elevation.png" width="600">

### Computations of elevations and variances:

This method calculates the satellite angle of elevation in the following stages:                     
Calculates the distance of receiver to the satellite (m) using pythagoras theorem.                   
Calculates the distance between the earth center and the satellite (m) using pythagoras theorem.     
Calculates the distance between the earth center and the receiver (m) using pythagoras theorem.      
These ranges make up a scalene triangle, where all ranges/sides are known.                                 
The low of cosines is used to calculate the angle about the receiver in degrees.                     
90 is subtracted from this angle to get the local elevation angle.                                   
The method then calculates the variance of the satellite at the calculated angle.                    
returns the variance as a float    

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Calculations/Elevations.png" width="600">

The code to achieve this is here:
```
    def elevation_variance_calculator(self) -> float:

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
        angle = radians((degrees(acos((ec_r ** 2 + r_s ** 2 - ec_s ** 2) / (2 * ec_r * r_s)))) - 90)

        return angle
    
    def variance(self):
        """""
        This method then calculates the variance of the satellite at the calculated angle.
        returns the variance as a float
        """""
        # Variance (uncertainty associated with the satellite) (m)
        variance = (self.l1_std ** 2) / (sin(self.elevation_variance_calculator()))
        
        return variance
        

```
Satellites ordered top to bottom: G24, G19, G18, G17, G15, G13, G12, G10

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Vectors/Degrees%20Vector.png" width="300"><img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Vectors/Radians%20Vector.png" width="300"><img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Vectors/Variances%20Vector.png" width="300">

(Elevation angles are representative) 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Diagrams/Satellite_Angles.png" width="600">

### l (Observations) vector 

This is the vector of raw observations.
UNITS: L1C (L1 cycles)

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Calculations/l.png" width="500">

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Vectors/Vector%20of%20observations%20(l)%20Vector.png" width="500">


### S (Single differencing) Matrix 

This matrix is used to generate a vector of single differences.

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Single%20Differencing%20(S)%20Matrix.png" width="500">

### Sl (Vector of single differences)

The dot product of the differencing matrix (S) and the vector of observations (l) generates the vector of single differences.

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Calculations/sl.png" width="500">

The following code was used compute this:
```
sl = S.dot(l)   
```
UNITS: L1C (L1 cycles)
<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Vectors/Vector%20of%20single%20differences%20(sl)%20Vector.png" width="500">

### D (Doube differencing) Matrix 

This matrix is used to generate values for the double differences of the observations matrix.

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Double%20Differencing%20(D)%20Matrix.png" width="500">

### DSl (Vector of Double Differences of Observations)

The dot product of the double differencing matrix (D) and the vector of single differences generates the vector of double differences of the observations

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Calculations/DSl.png" width="500">

The following code was used to compute this.

UNITS: L1C (L1 cycles)



Code:
```
Dsl = D.dot(sl)     
```

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Vectors/Double%20differences%20(Dsl)%20Vector.png" width="500">

### b (Observed - Computed)

The observed - computed (b) matrix was calculated in the following steps:

1. The phase ambiguity term (N) was subtracted from the double differenced vector of double differences (DSl):
In the code, this has been stored into the variable ``` DD_s_p_a``` This value is expressed in cycles

2. The range terms are computed from the satellite coordinates. These are then used to generate a value of the computed measurement. These are multiplied by 1/wavelength to convert from meters into cycles.

3. The observed measurement is subtracted from the computed measurement.

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Calculations/b.png" width="500">

UNITS: meters

The code to achieve this is here:
``` 
# This function is used to compute satelite - receiver ranges.
def distance(point_1: List[float], point_2: List[float]) -> float:
    return sqrt((point_2[0] - point_1[0])**2 +
                (point_2[1] - point_1[1])**2 +
                (point_2[2] - point_1[2])**2)
                
"""
wl - Wavelength 
brrs - Base receiver to reference satellite 
rrrs - Reference receiver to reference satellite
brcs - Base receiver to corresponding satellite 
rrcs - Reference receiver to corresponding satellite
DD_s_p_a - Vector of double differences, after phase ambiguity term subtracted.
"""    
brrs: floa = distance(ref_station, sat_ref)
rrrs = distance(rov_station, sat_ref)
brcs = distance(ref_station, corresponding_sat)
rrcs = distance(rov_station, corresponding_sat)
        
def calc_b_vector(wl: float, DD_s_p_a: float, brrs: float, rrrs: float, brcs: float, rrcs: float) -> float:                                        
    # observed - The vector of measured quantities                       
    o = DD_s_p_a                                                    
                                                                             
    # Computed                                                           
    c = (1 / wl) * (brrs - rrrs - brcs + rrcs)  
    return o - c                                                         
                                                                             
``` 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Vectors/Observed%20-%20Computed%20Vector.png" width="500">

### Cl (Observations covariance) Matrix

The covarince matrix of observations is calculated in the following steps:

Step 1: A 16 x 16 Identity matrix is generated.

Step 2: This matrix is multiplied by the vector of variances (Computed from satellite elevations)

This following code was used to make this computations: 
``` 
cl = np.eye(16, 16) * variance_vector  
``` 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Covariance%20matrix%20of%20observations%20(cl)%20Matrix.png" width="500">

### Cd (Covariance) Matrix 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Calculations/Cd.png" width="500">

``` 
def Cd_calculator(D, S, Cl):                                              
    return (((D.dot(S)).dot(Cl)).dot(transpose(S))).dot(transpose(D))      
``` 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/covariance%20matrix%20of%20the%20double%20differences%20(Cd)%20Matrix.png" width="500">

### Wd (Weight) Matrix 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Calculations/W.png" width="500">

The code to compute this is here:
``` 
Wd = linalg.inv(Cd)      
``` 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Weight%20(Wd)%20Matrix.png" width="500">

### A (Design) Martix

The design matrix was populated with the partial derivitives of the double difference observation equations with respect to the unknowns. 


<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Calculations/A.1.png" width="500">

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Calculations/A.2.png" width="500">


NOTE: Phase ambiguity terms are not included.

### The following code was used to perform these computations.

``` 
class DD:
    """"
    ref_station - Earth centric Cartesian Coordinates of the reference station [X, Y, Z]
    corresponding_sat - Earth centric cartesian Coordinates of the corresponding station [X, Y, Z]
    sat_ref - Satellite reference of cartesian Coordinates of the reference station [X, Y, Z]
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
        L1_f: float = 1575.42 * 1000000

        # Signal frequency of L2
        L2_f: float = 1227.6 * 1000000

        # Set to True by default
        if L1:
            wl = c / L1_f

        if L2:
            wl = c / L2_f

        # Error check the arguments
        assert len(ref_station) == len(rov_station) == len(corresponding_sat) == len(sat_ref)
        assert L1 != L2

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
                                                                           

``` 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Design%20(A)%20Matrix.png" width="500">

### Calculating the updates:

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Calculations/Xhat.png" width="500">

### ATWA Matrix 

``` 
def ATWA(A, W):                            
    return ((transpose(A)).dot(W)).dot(A)      
``` 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/ATWA%20Matrix.png" width="500">

### (ATWA)^-1 Matrix

``` 
linalg.inv(atwa)   
```

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/(ATWA)%5E-1%20Matrix.png" width="500">

## RESULT

### EPOCH 1
##### X update = 0.14130538m
##### Y update = -0.12730074m
##### Z update = 0.13206799m

##### Updated X Cooridnates: 4929605.54130538 
##### Updated Y Cooridnates: -29123.82730074 
##### Updated Z Cooridnates: 4033603.93206799 

### EPOCH 2
##### X update = 0.1410819m
##### Y update = -0.12767789m
##### Z update = 0.13099312m

##### Updated X Cooridnates: 4929605.5410819
##### Updated Y Cooridnates: -29123.82767789 
##### Updated Z Cooridnates: 4033603.93099312

### EPOCH 3
##### X update = 0.13911663m
##### Y update = -0.12739046m
##### Z update = 0.13271809m

##### Updated X Cooridnates: 4929605.53911663 
##### Updated Y Cooridnates: -29123.82739046
##### Updated Z Cooridnates: 4033603.93271809

## Position Solution as an average of 3 Epochs.

##### Averaged X Cooridnates: 4929605.540501303
##### Averaged Y Cooridnates: -29123.82745636333
##### Averaged Z Cooridnates: 4033603.9319264


## Quality Assessment:

The distance between the two pillars is approximatley 94.4m. 

Distance between Pillar 1A and Nominal coordinates for Pillar 3A: 
94.287m

Distance between Pillar 1A and updated coordinates for Pillar 3A
94.403m

This shows a 0.116 change in range.

#### Residuals:

Observed - Computed residuals are all within 99% confidence intervals for all three epochs. 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Graphs/Z%20Scores%20of%20Observed%20-%20Computed.png" width="500">





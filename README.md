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

Python 3

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

Satellite G24 has the highest elevation of 71 degrees. This satellite is used as the reference satellite.

It is important to initially calculate the elevation angles of each satelite. The error of a satellite is inversely proportial to elevation with respect to a local horizon. Low elevation satellites produce less reliable results, and this needs to be taken into account when formulating a weight matrix.


<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Diagrams/Satellite_Elevation.png" width="600">

(Elevation angles depicted here are to scale) 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Diagrams/Satellite_Angles.png" width="600">



## Matrices:

### l (Observations) vector 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Vector%20of%20observations%20(l)%20Matrix.png" width="500">

### A (Design) Matrix

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Design%20(A)%20Matrix.png" width="500">

### S (Single differencing) Matrix 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Single%20Differencing%20(S)%20Matrix.png" width="500">

### Sl (Vector of single differences)

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Vector%20of%20single%20differences%20(sl)%20Matrix.png" width="500">

### D (Doube differencing) Matrix 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Double%20Differencing%20(D)%20Matrix.png" width="500">

### DSl (Vector of Double Differences)

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Vector%20of%20Double%20differences%20(Dsl)%20Matrix.png" width="500">

### b (Observed - Computed)

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Vector%20of%20Double%20differences%20(Dsl)%20Matrix.png" width="500">

### Cl (Observations covariance) Matrix

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Covariance%20matrix%20of%20observations%20(cl)%20Matrix.png" width="500">

### Cd (Covariance) Matrix 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/covariance%20matrix%20of%20the%20double%20differences%20(Cd)%20Matrix.png" width="500">

### Wd (Weight) Matrix 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Weight%20(Wd)%20Matrix.png" width="500">

### A (Design) Martix

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/Design%20(A)%20Matrix.png" width="500">

### ATWA Matrix 

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/ATWA%20Matrix.png" width="500">

### (ATWA)^-1 Matrix

<img src="https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/Matrices/(ATWA)%5E-1%20Matrix.png" width="500">

### x^ result

### RESULT


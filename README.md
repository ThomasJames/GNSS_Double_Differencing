# GNSS Double Differencing

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

![Diagram aid](https://github.com/ThomasJames/GNSS_Double_Differencing/blob/master/DD_Diagram_aid.png)


## Prerequisites 

Python 3

The following libraries must be installed:

``` 
pip install numpy 
pip install matplotlib
pip install seaborn 
```

## Data

Original text file for data: ``` https://github.com/ThomasJames/GNSS_Data_(text).txt``` 
Data has been organised into .py file: ``` https://github.com/ThomasJames/GNSS_Double_Differencing/Data.py```   
NOTE: All data used in computations is derived from data from Data.py: Only data from the epoch ``` 2016 11 15 22 19 5```  is used.


## Method 

Satellite G24 has the highest elevation of 71 degrees azimuth. This satellite is used as the reference satellite.

#### Step 1:

Initial l vector was formed 

![](https://github.com/ThomasJames/Matrix_Images/I.png)

Single differences need to be derived:



Single differences were determined using a 8 x 16 differencing matrix D: 




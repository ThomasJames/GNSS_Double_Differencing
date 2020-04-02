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

Form a 'l' vector is formed using containing the 

A differencing matrix 's' is formed to calculate the single differences of 

![equation](<a href="https://www.codecogs.com/eqnedit.php?latex=l&space;=&space;[&space;\phi&space;^{G24}&space;_{1A}&space;\phi&space;^{G24}&space;_{3A}&space;\phi&space;^{G19}&space;_{1A}&space;\phi&space;^{G19}&space;_{3A}&space;\phi&space;^{G18}&space;_{1A}&space;\phi&space;^{G18}&space;_{3A}&space;\phi&space;^{G17}&space;_{1A}&space;\phi&space;^{G17}&space;_{3A}&space;\phi&space;^{G15}&space;_{1A}&space;\phi&space;^{G15}&space;_{3A}&space;\phi&space;^{G13}&space;_{1A}&space;\phi&space;^{G13}&space;_{3A}&space;\phi&space;^{G12}&space;_{1A}&space;\phi&space;^{G12}&space;_{3A}&space;\phi&space;^{G10}&space;_{1A}&space;\phi&space;^{G10}&space;_{3A}]&space;^T" target="_blank"><img src="https://latex.codecogs.com/gif.latex?l&space;=&space;[&space;\phi&space;^{G24}&space;_{1A}&space;\phi&space;^{G24}&space;_{3A}&space;\phi&space;^{G19}&space;_{1A}&space;\phi&space;^{G19}&space;_{3A}&space;\phi&space;^{G18}&space;_{1A}&space;\phi&space;^{G18}&space;_{3A}&space;\phi&space;^{G17}&space;_{1A}&space;\phi&space;^{G17}&space;_{3A}&space;\phi&space;^{G15}&space;_{1A}&space;\phi&space;^{G15}&space;_{3A}&space;\phi&space;^{G13}&space;_{1A}&space;\phi&space;^{G13}&space;_{3A}&space;\phi&space;^{G12}&space;_{1A}&space;\phi&space;^{G12}&space;_{3A}&space;\phi&space;^{G10}&space;_{1A}&space;\phi&space;^{G10}&space;_{3A}]&space;^T" title="l = [ \phi ^{G24} _{1A} \phi ^{G24} _{3A} \phi ^{G19} _{1A} \phi ^{G19} _{3A} \phi ^{G18} _{1A} \phi ^{G18} _{3A} \phi ^{G17} _{1A} \phi ^{G17} _{3A} \phi ^{G15} _{1A} \phi ^{G15} _{3A} \phi ^{G13} _{1A} \phi ^{G13} _{3A} \phi ^{G12} _{1A} \phi ^{G12} _{3A} \phi ^{G10} _{1A} \phi ^{G10} _{3A}] ^T" /></a>)



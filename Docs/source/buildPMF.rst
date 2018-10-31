.. highlight:: bash


buildPMF - Create Fortran 1D interpolator
*****************************************

Given a text file consisting of an array of states over a 1D set of points, create a fortran function that interpolates the states by computing the average of each state between two locations.  Typically, this is used to create a function for interpolating a 1D premixed flame solution from PREMIX or Cantera that can be linked with another code in 1-3 spatial dimensions for initializing data.


```
Usage:
```

Example:


# CalculateEpsilon
R code to calculate turbulent kinetic energy dissipation (epsilon) in running waters using Acoustic Doppler Velocimetry
Example data data ("ADV_ExampleFile") are also provided

code is written by Marcus Klaus, Swedish University of Agricultural Sciences (SLU), Umeå, Sweden
contact: marcus.klaus@posteo.net

last update: 2020-03-26
 
Please don't hesitate to contact me for any questions or if you find errors in the code!

Please note: There is no golden standard method for estimating epsilon in running waters
results are sensitive to many parameters, and must always be quality checked!

For use with different ADV, the code needs to be adjusted

The script consistes of three parts: 
(1) load directories / boundary data
(2) load functions
(3) process data

METHOD
Method is described in detail in 
Vingiani, F., Durighetto, N, Klaus, M., Schelker, J., Labasque, T. and Botter, W.: Evaluatingstream CO2 
outgassing via Drifting and Anchored chambers:results of a Flume experiment, for submission to Biogeosciences.
 
Code is based on methods in Zappa et al. 2003 and uses many improvements by Bluteau et al. 2011
Some code parts are adapted from the GDopp package (https://github.com/USGS-R/GDopp)

Bluteau, C., Jones, N., and Ivey, G.: Estimating turbulent kinetic energy dissipation using the inertial subrange method in environmental Flows, Limnology and Oceanography: Methods, 9, https://doi.org/10.4319/lom.2011.9.302, 2011.

Zappa, C. J., Raymond, P. A., Terray, E. A., and McGillis, W. R.: Variation in surface turbulence and the gas transfer velocity over a tidal cycle in a macro-tidal estuary, Estuaries, 26, 1401–1415, 2003.



# CalculateEpsilon
R code to calculate turbulent kinetic energy dissipation (epsilon) in running waters using Acoustic Doppler Velocimetry and the inertial subrange method.
Example data ("ADV_ExampleFile") are also provided

code is written by Marcus Klaus, Swedish University of Agricultural Sciences (SLU), Umeå, Sweden, 
contact: marcus.klaus@posteo.net

release date (code version 1_0): 2020-03-26

last update: (code version 1_1) 2020-11-11

The code consistes of four parts: 
(1) load directories / boundary data,
(2) load functions,
(3) process data,
(4) Quality check

# Please note
This code has been revised in response to reviewer comments; for details, see https://bg.copernicus.org/preprints/bg-2020-327/ 
Please use the latest version "CalcEpsilon1_1.R"

There is no golden standard method for estimating epsilon in running waters
results are sensitive to many parameters, and must always be quality checked!

This code is for a Nortek Vectrino 4 beam ADV. For use with different ADV systems, the code may need to be adjusted.

Please don't hesitate to contact me for any questions or if you find errors in the code!

# Method
The method is described in detail in

Vingiani, F., Durighetto, N, Klaus, M., Schelker, J., Labasque, T. and Botter, W.: Evaluatingstream CO2 
outgassing via Drifting and Anchored chambers:results of a Flume experiment. Biogeosciences Discussions https://bg.copernicus.org/preprints/bg-2020-327/ .
 
Thec ode is based on methods in Zappa et al. 2003 and uses many improvements by Bluteau et al. 2011. Some code parts are adapted from the GDopp package (https://github.com/USGS-R/GDopp) and the Matlab "Despike" functiom by Mori (2020).

Bluteau, C., Jones, N., and Ivey, G.: Estimating turbulent kinetic energy dissipation using the inertial subrange method in environmental Flows, Limnology and Oceanography: Methods, 9, https://doi.org/10.4319/lom.2011.9.302, 2011.

Mori, N. (2020). Despiking (https://www.mathworks.com/matlabcentral/fileexchange/15361-despiking), MATLAB Central File Exchange. Retrieved November 6, 2020.

Zappa, C. J., Raymond, P. A., Terray, E. A., and McGillis, W. R.: Variation in surface turbulence and the gas transfer velocity over a tidal cycle in a macro-tidal estuary, Estuaries, 26, 1401–1415, 2003.

For full reference list, see code and manuscript.



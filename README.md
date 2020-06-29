# CalculateEpsilon
R code to calculate turbulent kinetic energy dissipation in running waters using Acoustic Doppler Velocimetry
Example data data ("ADV_ExampleFile") are also provided

### GENERAL INFORMATION
### code is written by Marcus Klaus, Swedish University of Agricultural Sciences (SLU), Umeå, Sweden
### contact: marcus.klaus@posteo.net
### last update: 2020-03-26
### 
### Please don't hesitate to contact me for any questions or if you find errors in the code!
###
### Please note: There is no golden standard method for estimating epsilon in running waters
### results are sensitive to many parameters, and must always be quality checked!
### 
### For use with different ADV, the code needs to be adjusted
###
### The script consistes of three parts: 
### (1) load directories / boundary data
### (2) load functions
### (3) process data
###
###############################################################################################################
### METHOD
### Method is described in detail in 
### Vingiani, F., Durighetto, N, Klaus, M., Schelker, J., Labasque, T. and Botter, W.: Evaluatingstream CO2 
### outgassing via Drifting and Anchored chambers:results of a Flume experiment, for submission to Biogeosciences.
###  
### Code is based on methods in Zappa et al. 2003 and uses many improvements by Bluteau et al. 2011
### Some code parts are adapted from the GDopp package (https://github.com/USGS-R/GDopp)
###
### References:
###
### Bluteau, C., Jones, N., and Ivey, G.: Estimating turbulent kinetic energy dissipation using the inertial 
### subrange method in environmental Flows, Limnology and Oceanography: Methods, 9, https://doi.org/10.4319/lom.2011.9.302, 2011.
###
### Goring, D. and Nikora, V.: De-spiking Acoustic Doppler Velocimeter data., Journal of Hydraulic Engineering, 128, 117–126, 
### https://doi.org/10.1061/(ASCE)0733-9429(2002)128:1(117), 2002.
###
### Henjes, K., Taylor, P. K., and Yelland, M. J.: Effect of Pulse Averaging on Sonic Anemometer Spectra, 
### Journal of Atmospheric and Oceanic Technology, 16, 181–184, https://doi.org/10.1175/1520-0426(1999)016<0181:EOPAOS>2.0.CO;2, 1999.
###
### Kitaigorodskii, S. A. and Lumley, J. L.: Wave-Turbulence interactions in the Upper Ocean. Part I: 
### The Energy Balance of the Interacting Fields of Surface Wind Waves and Wind-Induced Three-Dimensional Turbulence, 
### Journal of Physical Oceanography, 13, 1977–1987, 
### https://doi.org/10.1175/1520-0485(1983)013<1977:WTIITU>2.0.CO;2,https://doi.org/10.1175/1520-0485(1983)013<1977:WTIITU>2.480 0.CO;2, 1983.
###
### Ruddick, B., Anis, A., and Thompson, K.: Maximum Likelihood Spectral Fitting: The Batchelor Spectrum, Journal of
### Atmospheric and Oceanic Technology - J ATMOS OCEAN TECHNOL, 17, 1541–1555, https://doi.org/10.1175/15200426(2000)017<1541:MLSFTB>2.0.CO;2, 2000.
###
### Zappa, C. J., Raymond, P. A., Terray, E. A., and McGillis, W. R.: Variation in surface turbulence and the 
### gas transfer velocity over a tidal525 cycle in a macro-tidal estuary, Estuaries, 26, 1401–1415, 2003.
###


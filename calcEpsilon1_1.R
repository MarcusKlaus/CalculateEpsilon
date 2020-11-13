###############################################################################################################
### This R code calculates turbulent kinetic energy dissipation rate (epsilon) from 4D Nortek Vectrino ADV data
###############################################################################################################
### code written by Marcus Klaus, Postdoc, Swedish University of Agricultural Sciences (SLU), Umeå, Sweden
### contact: marcus.klaus@posteo.net
### last update: 2020-11-11
### 
### Please don't hesitate to contact me for any questions or if you find errors in the Code!
###
### Please note: There is no golden standard method for estimating epsilon in rivers
### results are sensitive to many parameters, and must always be quality checked!
###
### I also provide example data files ("ADV_ExampleFile") collected using a 4D Nortek Vectrino ADV 
### in a flume at Lunz Mesocosm facility, Lunz am See, Austria
### 
### Please note: For use with different ADV, the code needs to be adjusted
###
### The script consistes of four parts: 
### (1) load directories / boundary conditions
### (2) load functions
### (3) process data
### (4) Quality check of results
###
###############################################################################################################
### METHOD
### Method is described in detail in
###  
### Vingiani, F., Durighetto, N, Klaus, M., Schelker, J., Labasque, T. and Botter, W.: Evaluatingstream CO2 
### outgassing via Drifting and Anchored chambers:results of a flume experiment, Biogeosciences Discussions. https://bg.copernicus.org/preprints/bg-2020-327/
###  
###
### The code is based on methods in Zappa et al. 2003 and uses many improvements by Bluteau et al. 2011
### Some code parts are adapted from the GDopp R package (https://github.com/USGS-R/GDopp) and the Matlab "Despike" function by Mori (2020)
###
### References:
###
### Baker, M.A. and Gibson C.H.: Sampling Turbulence in the Stratified Ocean: Statistical Consequences of Strong Intermittency. Journal of Physical Oceanography 17, 1817-1836. 
### https://doi.org/10.1175/1520-0485(1987)017<1817:STITSO>2.0.CO;2, 1987
### 
### Bluteau, C., Jones, N., and Ivey, G.: Estimating turbulent kinetic energy dissipation using the inertial 
### subrange method in environmental ?ows, Limnology and Oceanography: Methods, 9, https://doi.org/10.4319/lom.2011.9.302, 2011.
###
### Doroudian, B, Bagherimiyab, F. and Lemmin, U.: Improving the accuracy of four-receiver acoustic Doppler velocimeter (ADV) measurements in turbulent boundary layer flow
### Limnol. Oceanogr.: Methods 8, 575–591. https://doi.org/10.4319/lom.2010.8.0575, 2010.
### 
### Thomas Foken (2005): Micrometeorology. Springer.
###
### Goring, D. and Nikora, V.: De-spiking Acoustic Doppler Velocimeter data., Journal of Hydraulic Engineering, 128, 117–126, 
### https://doi.org/10.1061/(ASCE)0733-9429(2002)128:1(117), 2002.
###
### Henjes, K., Taylor, P. K., and Yelland, M. J.: Effect of Pulse Averaging on Sonic Anemometer Spectra, 
### Journal of Atmospheric and Oceanic Technology, 16, 181–184, https://doi.org/10.1175/1520-0426(1999)016<0181:EOPAOS>2.0.CO;2, 1999.
###
### Kitaigorodskii, S. A. and Lumley, J. L.: Wave-Turbulence interactions in the Upper Ocean. Part I: 
### The Energy Balance of the Interacting Fields of Surface Wind Waves and Wind-Induced Three-Dimensional Turbulence, 
### Journal of Physical Oceanography, 13, 1977–1987, https://doi.org/10.1175/1520-0485(1983)013<1977:WTIITU>2.0.CO;2,https://doi.org/10.1175/1520-0485(1983)013<1977:WTIITU>2.480 0.CO;2, 1983.
###
### Lohrmann, Atle, Ramon Cabrera, and Nicholas C. Kraus: Acoustic-Doppler velocimeter (ADV) for laboratory use. 
### Fundamentals and advancements in hydraulic measurements and experimentation, pp. 351-365. ASCE, 1994.
###
### Mori, N. (2020). Despiking (https://www.mathworks.com/matlabcentral/fileexchange/15361-despiking), MATLAB Central File Exchange. Retrieved November 6, 2020.
###
### Mori, N., T. Suzuki and S. Kakuno: Noise of acoustic Doppler velocimeter data in bubbly flow, Journal of Engineering Mechanics, American Society of Civil Engineers, Volume 133, Issue 1, pp.122-125, 
### https://doi.org/10.1061/(ASCE)0733-9399(2007)133:1(122), 2007
###
### Ruddick, B., Anis, A., and Thompson, K.: Maximum Likelihood Spectral Fitting: The Batchelor Spectrum, Journal of
### Atmospheric and Oceanic Technology - J ATMOS OCEAN TECHNOL, 17, 1541–1555, https://doi.org/10.1175/15200426(2000)017<1541:MLSFTB>2.0.CO;2, 2000.
###
### Zappa, C. J., Raymond, P. A., Terray, E. A., and McGillis, W. R.: Variation in surface turbulence and the 
### gas transfer velocity over a tidal cycle in a macro-tidal estuary, Estuaries, 26, 1401–1415, 2003.
###
###############################################################################################################
### INPUTS: 
# ADV files (.dat and .hdr)
# ADV geometry (roll, pitch and heading of sensor)
# TransmitDepth		distance between water surface and ADV transmit transducer (m)
# WaterWidth		stream width (m)
# heading			ADV probe heading (degrees)
# pitch			ADV probe pitch (degrees)
# roll			ADV probe roll (degrees)
#
### OUTPUTS: 

# (1) plots with velocity time series, and spikes highlighted
# (2) plot with velocity spectra
# (3) output.csv-file which contains the following columns:
# Variable name:		Explanation:
# filename			Filename
# DataTime			Date and time
# VelocityRange		Vectrino setting of Velocity range (m/s)
# TransmitLength		Vectrino setting of Transmit length (mm)
# SamplingLength		Vectrino setting of Sampling volume length scale (mm)
# Distance			Distance between sensor and "wall"
# Temperature		Water temperature (dC)
# VelocityScaling		Vectrino setting of Velocity scaling
# CoordinateSystem	Coordinate system of data aquisition
# WaterDepth		Water depth (m)
# SensorDepth		Distance of ADV sampling volume from water surface (m)
# TKE				Turbulent kinetic energy (m2/s2)
# Re				Reynolds number
# Fr				Froude number
# eD				Energy dissipation (by drag) (m2/s3)
# eS				Energy dissipation (by shear) (m2/s3)
# V				mean flow velocity in u direction (m/s)
# signal.rat.X		signal to noise ratio for X direction
# signal.rat.Y		signal to noise ratio for Y direction
# signal.rat.Z1		signal to noise ratio for Z1 direction
# signal.rat.Z2		signal to noise ratio for Z2 direction
# correlation.X		correlation for X direction
# correlation.Y		correlation for Y direction
# correlation.Z1		correlation for Z1 direction
# correlation.Z2		correlation for Z2 direction
# FrozenTurb.u		Metric to test Frozen turbulence hypothesis (u direction)
# FrozenTurb.v		Metric to test Frozen turbulence hypothesis (v direction)
# FrozenTurb.w1		Metric to test Frozen turbulence hypothesis (w1 direction)
# FrozenTurb.w2		Metric to test Frozen turbulence hypothesis (w2 direction)
# propDespiked.u		Proportion of data removed by Despiking (u direction)
# propDespiked.v		Proportion of data removed by Despiking (v direction)
# propDespiked.w1		Proportion of data removed by Despiking (w1 direction)
# propDespiked.w2		Proportion of data removed by Despiking (w2 direction)
# sdVel.u			Standard deviation of velocity (u direction) (m/s)
# sdVel.v			Standard deviation of velocity (v direction) (m/s)
# sdVel.w1			Standard deviation of velocity (w1 direction) (m/s)
# sdVel.w2			Standard deviation of velocity (w2 direction) (m/s)
# epsilonLS:u		Least square estimate of TKE dissipation rate (m2/s3) (u direction); based on unfiltered spectrum
# epsilonLS:v		Least square estimate of TKE dissipation rate (m2/s3) (v direction); based on unfiltered spectrum
# epsilonLS:w1		Least square estimate of TKE dissipation rate (m2/s3) (w1 direction); based on unfiltered spectrum
# epsilonLS:w2		Least square estimate of TKE dissipation rate (m2/s3) (w2 direction); based on unfiltered spectrum
# epsilonMLE:u		Maximum Likelihood estimate of TKE dissipation rate (m2/s3) (u direction), mean; based on unfiltered spectrum
# epsilonMLE:v		Maximum Likelihood estimate of TKE dissipation rate (m2/s3) (v direction), mean; based on unfiltered spectrum
# epsilonMLE:w1		Maximum Likelihood estimate of TKE dissipation rate (m2/s3) (w1 direction), mean; based on unfiltered spectrum
# epsilonMLE:w2		Maximum Likelihood estimate of TKE dissipation rate (m2/s3) (w2 direction), mean; based on unfiltered spectrum
# epsilonMLElwr:u		Maximum Likelihood estimate of TKE dissipation rate (m2/s3) (u direction), lower boundary of 95% confidence interval; based on unfiltered spectrum
# epsilonMLElwr:v		Maximum Likelihood estimate of TKE dissipation rate (m2/s3) (v direction),  lower boundary of 95% confidence interval; based on unfiltered spectrum
# epsilonMLElwr:w1	Maximum Likelihood estimate of TKE dissipation rate (m2/s3) (w1 direction),  lower boundary of 95% confidence interval; based on unfiltered spectrum
# epsilonMLElwr:w2	Maximum Likelihood estimate of TKE dissipation rate (m2/s3) (w2 direction),  lower boundary of 95% confidence interval; based on unfiltered spectrum
# epsilonMLEupr:u		Maximum Likelihood estimate of TKE dissipation rate (m2/s3) (u direction), upper boundary of 95% confidence interval; based on unfiltered spectrum
# epsilonMLEupr:v		Maximum Likelihood estimate of TKE dissipation rate (m2/s3) (v direction),  upper boundary of 95% confidence interval; based on unfiltered spectrum
# epsilonMLEupr:w1	Maximum Likelihood estimate of TKE dissipation rate (m2/s3) (w1 direction),  upper boundary of 95% confidence interval; based on unfiltered spectrum
# epsilonMLEupr:w2	Maximum Likelihood estimate of TKE dissipation rate (m2/s3) (w2 direction),  upper boundary of 95% confidence interval; based on unfiltered spectrum
# realr2:u			Coefficient of determination of maximum likelihood predicted vs. observed spectrum (u direction), based on unfiltered spectrum
# realr2:v			Coefficient of determination of maximum likelihood predicted vs. observed spectrum (v direction), based on unfiltered spectrum
# realr2:w1			Coefficient of determination of maximum likelihood predicted vs. observed spectrum (w1 direction), based on unfiltered spectrum
# realr2:w2			Coefficient of determination of maximum likelihood predicted vs. observed spectrum (w2 direction), based on unfiltered spectrum
# MAD:u			Maximum absolute deviation of maximum likelihood predicted vs. observed spectrum (u direction), based on unfiltered spectrum
# MAD:v			Maximum absolute deviation of maximum likelihood predicted vs. observed spectrum (v direction), based on unfiltered spectrum
# MAD:w1			Maximum absolute deviation of maximum likelihood predicted vs. observed spectrum (w1 direction), based on unfiltered spectrum
# MAD:w2			Maximum absolute deviation of maximum likelihood predicted vs. observed spectrum (w2 direction), based on unfiltered spectrum
# lower.k:u			selected lower wave number threshold (rad/m) for inertial subrange (u direction), based on unfiltered spectrum
# lower.k:v			selected lower wave number threshold (rad/m) for inertial subrange (v direction), based on unfiltered spectrum
# lower.k:w1		selected lower wave number threshold (rad/m) for inertial subrange (w1 direction), based on unfiltered spectrum
# lower.k:w2		selected lower wave number threshold (rad/m) for inertial subrange (w2 direction), based on unfiltered spectrum
# upper.k:u			selected upper wave number threshold (rad/m) for inertial subrange (u direction), based on unfiltered spectrum
# upper.k:v			selected upper wave number threshold (rad/m) for inertial subrange (v direction), based on unfiltered spectrum
# upper.k:w1		selected upper wave number threshold (rad/m) for inertial subrange (w1 direction), based on unfiltered spectrum
# upper.k:w2		selected upper wave number threshold (rad/m) for inertial subrange (w2 direction), based on unfiltered spectrum
# f:u				frequency associated with kolmogorov scale (above which spectrum shows noise) (u direction), based on unfiltered spectrum
# f:v				frequency associated with kolmogorov scale (above which spectrum shows noise) (v direction), based on unfiltered spectrum
# f:w1			frequency associated with kolmogorov scale (above which spectrum shows noise) (w1 direction), based on unfiltered spectrum
# f:w2			frequency associated with kolmogorov scale (above which spectrum shows noise) (w2 direction), based on unfiltered spectrum

###############################################################################################################

################################################################################################
################################################################################################
################################################################################################
## Part 1: setting directories and boundary conditions
################################################################################################
################################################################################################
################################################################################################

###############################################
### set boundary conditions not measured by ADV
###############################################

## use conditions for ADV example file 
TransmitDepth = 0.02
WaterWidth = 0.4

## define geometrical position of ADV
position_data <- NULL
## heading
position_data$heading <- 0
## pitch
position_data$pitch <- 56
## roll
position_data$roll <- 45

## flume channel slope
slope=0.0025

## should velocity time series and spectra be plotted?
plotting = TRUE

## number of bins over which spectrum is averaged. 
## this is used for identifying bounds of inertial subrange to make sure it spans at least one decade in spectral density
## choosing bin size of 50 to remove noise and retain meaningful patterns
SpAve=51


###############################################
## set directory to folder with ADV files
###############################################
 
setwd("C:/MARCUS/WORK/Applications/Aquacosm_TA/Data/GitHUB_calcEpsilon/revision/")
folder.nm <- getwd()

################################################################################################
################################################################################################
################################################################################################
## Part 2: Loading functions
################################################################################################
################################################################################################
################################################################################################

###############################################
### load packages
###############################################

library(oce)
library(forecast)
library(signal)
library(emdbook)
library(gdata)
library(pracma)
library(lattice)

###############################################
### load functions for hydrodynamic indices
###############################################

## turbulent kinetic energy dissipation [m2/s3] (Raymond et al. 2012)
#(driven by drag)
Edf <- function(S,V){
Ed=9.81*S*V
return(Ed)
}
#(driven by shear) ## NOTE: this is only valid for smooth flow with small roughness!
Esf <- function(S,D,W){
Rh=D*W/(W+2*D)
us=sqrt(9.81*Rh*S)
Es=us^3/D
return(Es)
}

## Froude number
Frf <- function(V,D){
Fr=V/sqrt(9.81*D)
return(Fr)
}

## Reynolds number
Ref <- function(V,D,kinvis){
R=D*0.4/(2*D+0.4) # hydraulic radius
Re=V*R/kinvis
return(Re)
}

## Turbulent kinetic energy (TKE)
TKEf <- function(su,sv,sw){
TKE=1/2*sqrt(su^2+sv^2+sw^2)
return(TKE)
}


## geometric standard deviation
geomSD <- function(x,geoMean){
exp( sqrt( sum((log(x/geoMean))^2)/(length(x)-1) ) )
}


##############################################################################################
### Function for loading ADV data (adapted from https://github.com/USGS-R/GDopp )
##############################################################################################

load_adv <- function(file.nm, folder.nm){
  
  drop.cols <- c('checksum')
  adv.dat.names <- c('burst.num','ensemble.num','velocity.X','velocity.Y','velocity.Z1','velocity.Z2',
                     'amplitude.X','amplitude.Y','amplitude.Z1','amplitude.Z2',
                     'signal.rat.X','signal.rat.Y','signal.rat.Z1','signal.rat.Z2',
                     'correlation.X','correlation.Y','correlation.Z1','correlation.Z2')
  file.loc <- file.path(folder.nm,file.nm)
  data.adv <- read.table(file.loc)
  
  names(data.adv) <- adv.dat.names
  
  rmv.i <- data.adv$checksum==1
  if (any(rmv.i)){
    warning('some checksum failures. Removing errant values')
    data.adv <- data.adv[!rmv.i, ]
  }
  
  data.adv <- data.adv[,!(names(data.adv) %in% drop.cols)]
  
  return(data.adv)
}


##############################################################################################
### Function to derive epsilon from ADV flow velocity time series
############################################################################################## 
# Solves for epsilon from ADV data in the inertial subrange following Zappa et al. 2003, Bluteau et al. 2011
# Requires ADV dataset loaded by function "load_adv", and further instrumental and environmental conditions
# Optional: - Time series filtering (filt=TRUE)
#           - plot spectra (diagnostic=TRUE)
# Returns summary table with epsilon, different quality measures and other useful hydrodynamic parameters 


fit_epsilon <- function(dataset,freq=200, WaterDepth=WaterDepth, SensorDepth=SensorDepth,SamplingLength=SamplingLength,
diagnostic = FALSE,alpha=1.5*18/55*c(1,1.33,1.33,1.33),filt=TRUE,n=12,fmax=50,SpAve=11){
  
  ## select flow components
  u <- dataset[,1]
  v <- dataset[,2]
  w1 <- dataset[,3]
  w2 <- dataset[,4]

  ## low-pass filtering to remove high-frequency noise
  if(filt==TRUE){ 
  ## if fmax is higher than nyquist f than set filter to 0.95
  for(k in 1:4){
  if(fmax[k]>95){fmax[k]=95}
  }
  ## frequency expressed in terms of nyquist frequency and [0,1]
  uf <- filter(butter(n, fmax[1]/(0.5*freq), type="low"), u)
  vf <- filter(butter(n, fmax[2]/(0.5*freq), type="low"), v)
  w1f <- filter(butter(n, fmax[3]/(0.5*freq), type="low"), w1)
  w2f <- filter(butter(n, fmax[4]/(0.5*freq), type="low"), w2)
  }else{
  uf <- u
  vf <- v
  w1f <- w1
  w2f <- w2
  }

  ## make time series
  uts <- ts(uf, frequency=freq)
  vts <- ts(vf, frequency=freq)
  w1ts <- ts(w1f, frequency=freq)
  w2ts <- ts(w2f, frequency=freq)
 
  ## make Welch spectra using default 8 segments with 50% overlap and tapered with a Hamming window
  wu <- pwelch(uts, plot=FALSE)
  wv <- pwelch(vts, plot=FALSE)
  ww1 <- pwelch(w1ts, plot=FALSE)
  ww2 <- pwelch(w2ts, plot=FALSE) 
  wavenum.spectra.u <- wu$spec
  wavenum.spectra.v <- wv$spec
  wavenum.spectra.w1 <- ww1$spec
  wavenum.spectra.w2 <- ww2$spec 

  ##########################################
  ## mean advection velocity; based on vector rotated (u,v,w)->(u',0,0)
  ## calculated as in Bluteau et al. 2011, using "project" function  
  v.mn <- velocity_calc(project(dataset))

  ##########################################  
  ## calculate characteristic length scales

  # conversion from frequency to wavenumber space (rad/m)
  wavenum <- 2*pi*wu$freq/v.mn 
  # length scale of sampling volume; divide by 1000 to convert mm to m
  vol.k = 2*pi*(1/SamplingLength*1000) 
  # wave number corresponding to stream width
  width.k = 2*pi*(1/WaterWidth)
  # wave number corresponding to length scale of the largest isotropic eddies 
  # assuming this is the distance to the water surface
  depth.k = 2*pi*(1/SensorDepth)
  # noise limit based on fmax (Hz), also used for time series filtering 
  noise.k = 2*pi*fmax/v.mn 

  ## pulse average spectrum following Henjes et al. 1999
  pulse.averaging <- function(w.s,f,del.t){  
    denom1 <- (sin(pi*f*del.t)/(pi*f*del.t))^2 
    denom2 <- (f/(max(f)*2-f))^(5/3)*((sin(pi*(max(f)*2-f)*del.t))/(pi*(max(f)*2-f)*del.t))^2 
    pc.w = w.s/(denom1+denom2)
    return(pc.w)  
  }

  ## pulse averaging spectra
  w.s <- NULL
  w.s[[1]] <- pulse.averaging(wavenum.spectra.u,f=wu$freq,del.t=1/freq)
  w.s[[2]] <- pulse.averaging(wavenum.spectra.v,f=wv$freq,del.t=1/freq)
  w.s[[3]] <- pulse.averaging(wavenum.spectra.w1,f=ww1$freq,del.t=1/freq)
  w.s[[4]] <- pulse.averaging(wavenum.spectra.w2,f=ww2$freq,del.t=1/freq)

####################################################################################   
## find inertial subrange within which epsilon is calculated
## lower wavenumber limit is defined by stream depth, upper wavenumber limit is defined by ADV sampling volume
## fit spectrum for a series of candidate windows which where spectral density falls of by at least half a decade  
## loop through candidate windows to find window with best model fit, i.e. lowest mean average deviation (MAD) (loosely following Bluteau et al. 2011)

## make output vectors
lower.k.s <- matrix(NA,1000,4)
upper.k.s <- matrix(NA,1000,4)
ratio.k.s <- matrix(NA,1000,4)
realr2.s <- matrix(NA,1000,4)
MAD.s <- matrix(NA,1000,4)
epsilonMLEmean.s <- matrix(NA,1000,4)
epsilonMLElwr.s <- matrix(NA,1000,4)
epsilonMLEupr.s <- matrix(NA,1000,4)
epsilonMLE <- NULL
epsilonMLElwr <- NULL
epsilonMLEupr <- NULL
realr2 <- NULL
MAD <- NULL
lower.k <- NULL
upper.k <- NULL
ratio.k <- NULL

## loop through all h ADV beam components
for (h in 1:4){

# set start values
s=2
# set upper window limit to wave number associated with ADV sampling volume
upper.k.s[1,h] <- vol.k
lower.k.s[1,h] <- depth.k

# loop through candidate windows
while (upper.k.s[s-1,h] <= vol.k){
lower.k.s[s,h] <- wavenum[ which(wavenum >= depth.k) ][s]

### find upper bound; there must be minimum a decade in spectral density between upper and lower bound (Bluteau et al. 2011); 
## average spectral density over "spAve" bins to reduce noise
upper.k.s[s,h] <- wavenum[min(which(runmed(w.s[[h]],SpAve)[which(lower.k.s[s,h] == wavenum)]/runmed(w.s[[h]],SpAve) > 10  ))] 

## take note of the ratio  of spectral density between upper and lower bound
if(is.na(upper.k.s[s,h])==FALSE){
ratio.k.s[s,h] <- runmed(w.s[[h]],SpAve)[which(wavenum==lower.k.s[s,h])]/runmed(w.s[[h]],SpAve)[which(wavenum==upper.k.s[s,h])]
}

## if there is no upper boundary found, stop loop
if(is.na(upper.k.s[s,h])==TRUE){
break
}
if(upper.k.s[s,h]> vol.k){
break
}
if(lower.k.s[s,h] >= upper.k.s[s,h]){
break
}

if(s == 1000){
break
}

use.i.s <- lower.k.s[s,h] <= wavenum & wavenum <= upper.k.s[s,h] 

# select spectral data in candidate window
obsK <- w.s[[h]][use.i.s]

#########################################################
### Fit epsilon using MLE (Bluteau et al. 2011, eq. 15)
n=length(obsK)
d=wu$df
rm(a,fa)
logL <- NULL

## generate prior for potential epsilon values
ep <- lseq(0.0000001,1,1000)

## loop through this range
for(e in 1:length(ep)){
predK <- alpha[h]*(ep[e]^(2/3))*(wavenum[use.i.s])^(-5/3)
a=d*obsK/predK

# probability density function of the chi-squared distributed sample
fa <- (a^((d-2)/2)*exp(-a/2))/(2^(d/2)*gamma(d/2))
logL[e] = n*log(d)-sum(log(predK))+sum(log(fa))
} # end e loop

epsilonMLEmean.s[s,h] <- ep[which(logL==max(logL))]

# estimate variance of epsilon, i.e. the curvature of logL at its maximum (Bluteau et al. 2011, eq. 16)
if(which(logL==max(logL))<3){
vare <- -1/(diff(diff(logL))/diff(diff(ep)))[1]
}else{
vare <- -1/(diff(diff(logL))/diff(diff(ep)))[which(logL==max(logL))-2] ## -2 to account for difference calculations
}

# estimate asymptotic 95% confidence interval of epsilon
epsilonMLElwr.s[s,h] <- ep[which(logL==max(logL))] - 1.96*sqrt(vare)
epsilonMLEupr.s[s,h] <- ep[which(logL==max(logL))] + 1.96*sqrt(vare)


### fit data using MLE and evaluate fit  (Bluteau et al. 2011, Ruddick et al. 2000)

    ## predict spectral data
    predK <- alpha[h]*(epsilonMLEmean.s[s,h]^(2/3))*wavenum[use.i.s]^(-5/3)

    ## calculate R2 of model fit 
    SStot <- sum((obsK-mean(obsK))^2)
    SSres <- sum((obsK-predK)^2)
    #SSreg <- sum((predK-mean(obsK))^2)
    realr2.s[s,h] <- 1-(SSres/SStot)

    ## calulate maximum absolute deviation (MAD)
    MAD.s[s,h] <- 1/length(obsK)* sum( abs(obsK/predK-mean(obsK/predK)) ) 
    ## reject MAD?
    if(is.na(MAD.s[s,h])==FALSE){
    if(MAD.s[s,h] > 2*(2/d)^(1/2)){
    MAD.s[s,h] <- NA
    }
    }

s=s+1
} ## close s loop

## select candidate window for which MAD is minimal (Bluteau et al. 2011)
if(s>3){
id <- which(MAD.s[1:(s-2),h]==min(MAD.s[1:(s-2),h],na.rm=TRUE))
}else{
id = 2
} ## close s loop

## set epsilon estimates to NA if there is no valid estimate (MAD did not pass quality criteria)
if(length(id)==0){
epsilonMLE[h] <- NA
epsilonMLElwr[h] <- NA
epsilonMLEupr[h] <- NA
realr2[h] <- NA
}else{
## select final (best) epsilon estimate and save estimate, goodness of fit measures and confidence interval
epsilonMLE[h] <- epsilonMLEmean.s[id,h]
epsilonMLElwr[h] <- epsilonMLElwr.s[id,h]
epsilonMLEupr[h] <- epsilonMLEupr.s[id,h]
realr2[h] <- realr2.s[id,h]
MAD[h] <- MAD.s[id,h]
lower.k[h] <- lower.k.s[id,h]
upper.k[h] <- upper.k.s[id,h]
ratio.k[h] <- ratio.k.s[id,h]
}

} ## close h loop

### kolmogorov length scale (threshold for noise)
k.vis <- get_kin_viscosity(Temperature)
nk <- (k.vis^3/epsilonMLE)^(1/4)
nk.k <- 1/(10*nk)
	
## plot spectra for all ENU directions

  if (diagnostic){
    
    dev.new(width=8, height=8)
    layout(matrix(c(1:1), 1, 1, byrow = TRUE))
    par(mar=c(4.5,4.5,0.5,1),las=1,cex=1.3,mgp=c(3.1, 0.8, 0))

    plot(wavenum,w.s[[1]],log='xy',type="l",xlab=expression(paste("wave number [rad ", m^-1,"]")),ylab=expression(paste("Power spectral density [( ",m^2, s^-2,")","(rad ",m^-1,")",.^-1,"]")),ylim=c(10^-8,10^-2),xaxt="n",yaxt="n")
    
    axis(side=2,c(seq(0.00000001,0.0000001,0.00000001),seq(0.0000001,0.000001,0.0000001),seq(0.000001,0.00001,0.000001),seq(0.00001,0.0001,0.00001),seq(0.0001,0.001,0.0001),seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,1,0.1),seq(1,10,1),seq(10,100,10),seq(100,1000,100),seq(1000,10000,1000)),labels=FALSE)
    axis(side=2,c(0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000),lwd=2)
    axis(side=1,c(seq(0.00001,0.0001,0.00001),seq(0.0001,0.001,0.0001),seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,1,0.1),seq(1,10,1),seq(10,100,10),seq(100,1000,100),seq(1000,10000,1000)),labels=FALSE)
    axis(side=1,c(0.01,0.1,1,10,100,1000,10000),labels=c(0.01,0.1,1,10,100,1000,10000),lwd=2)

    points(wavenum,w.s[[2]],col='red',cex=.2,type="l")
    points(wavenum,w.s[[3]],col='green',cex=.2,type="l")
    points(wavenum,w.s[[4]],col='blue',cex=.2,type="l")
    #abline(v=width.k,lty=2)
    #text(width.k*1.15,5*10^-3,"width",srt=90)
    abline(v=depth.k,lty=2)
    text(depth.k*1.15,5*10^-3,"length",srt=90)
    abline(v=vol.k,lty=3)
    text(vol.k*1.15,5*10^-3,"vol",srt=90)


    ## add lines based on which epsilon is estimated
    for (p in 1:4){
    ip <- wavenum>lower.k[p] & wavenum<upper.k[p]
    lines(wavenum[ip],alpha[p]*(epsilonMLE[p]^(2/3))*wavenum[ip]^(-5/3),col=p+4,lwd=2)
    # add confidence intervals
    #lines(wavenum[ip],alpha[p]*(epsilonMLElwr[p]^(2/3))*wavenum[ip]^(-5/3),col=p+4,lwd=1)
    #lines(wavenum[ip],alpha[p]*(epsilonMLEupr[p]^(2/3))*wavenum[ip]^(-5/3),col=p+4,lwd=1)
    }

    text(min(wavenum),3*10^-7,paste0("R2u=",round(realr2[1],2)))
    text(min(wavenum),10^-7,paste0("R2v=",round(realr2[2],2)))
    text(min(wavenum),3*10^-8,paste0("R2w1=",round(realr2[3],2)))
    text(min(wavenum),10^-8,paste0("R2w2=",round(realr2[4],2)))

    text(min(wavenum),3*10^-5,paste0("eu=",round(epsilonMLE[1],4)))
    text(min(wavenum),10^-5,paste0("ev=",round(epsilonMLE[2],4)))
    text(min(wavenum),3*10^-6,paste0("ew1=",round(epsilonMLE[3],4)))
    text(min(wavenum),10^-6,paste0("ew2=",round(epsilonMLE[4],4)))


### kolmogorov length scale (threshold for noise)
abline(v=nk.k,lty=3,col=c(1:4))

legend("topright",c("u","v","w1","w2"),col=c(1,2,3,4),lwd=1,bg="white")

}

## frequency associated with kolmogorov scale; here I use it as input for filtering the spectrum
f= nk.k*v.mn/(2*pi)

## generate output
out <- rbind(epsilonMLE,epsilonMLElwr,epsilonMLEupr,realr2,MAD,lower.k,upper.k,f,ratio.k)
colnames(out) <- c("u","v","w1","w2")
return(out)
}


##########################################
# get kinematic viscosity
##########################################

# code from GDopp package ( https://github.com/USGS-R/GDopp )

# temperature in degC of water

get_kin_viscosity <- function(temperature=20) {
  # from Mays 2005, Water Resources Engineering
  tempTable <- seq(0,100,by=5)
  # table in m2/s E-6
  visTable <- c(1.792,1.519,1.308,1.141,1.007,0.897,
                0.804,0.727,0.661,0.605,0.556,0.513,0.477,0.444,
                0.415,0.39,0.367,0.347,0.328,0.311,0.296)
  v <- data.frame(approx(tempTable,visTable,xout = temperature))[2]
  v <- v*1e-6
  return(v$y)
}


##########################################
# coordinate transformation
##########################################

# code adapted from GDopp package ( https://github.com/USGS-R/GDopp )
# original code for coordinate transformation is based on 3*3 matrix (from http://www.nortek-
# as.com/en/knowledge-center)
# Marcus Klaus received 4*4 tilt and heading matrices by personal communication with Nortek on 2020-01-17

# transform velocity data between beam coordinates (along ADV beam direction) and ENU coordinates, where
# E = East-West component, N = North-South component, U = Up-Down component

# transformation follows Lohrman et al. (1994)


coord_transform <- function(trans_matrix, data_v, position_data,Z="Z1"){
    
  trans_matrix1 <- trans_matrix

  x <- data_v$velocity.X
  y <- data_v$velocity.Y
  z1 <- data_v$velocity.Z1
  z2 <- data_v$velocity.Z2

  xyz <- matrix(data=c(x,y,z1,z2), ncol = 4) 

  heading <- position_data$heading
  pitch <- position_data$pitch
  roll <- position_data$roll 
  res_trans <- resultant_trans(trans_matrix1, heading, pitch, roll)

  ENU <- apply(X = t(xyz), MARGIN = 2, FUN = xyz_2_enu, 
               trans_matrix = trans_matrix1, res_trans = res_trans)#orig

  ENU.df <- data.frame(t(ENU))
  ## change order
  ENU.df <- ENU.df[,c(2,1,3,4)]
  names(ENU.df) <- c('North','East','Up1','Up2')
  return(ENU.df)
}

resultant_trans <- function(trans_matrix1, heading, pitch, roll){

  hh <- pi * (heading - 90) / 180
  pp <- pi * pitch / 180
  rr <- pi * roll / 180
  
  # Make heading matrix
  H <- matrix(data = c(cos(hh), sin(hh), 0, 0, 
				-sin(hh), cos(hh), 0, 0,  
				0, 0, 1, 0,  
				0, 0, 0, 1), ncol = 4)

  # Make tilt matrix
  p_data <- c(cos(pp), -sin(pp)*sin(rr), -cos(rr)*sin(pp), -cos(rr)*sin(pp),
		0, cos(rr), -sin(rr), -sin(rr),
           sin(pp), sin(rr)*cos(pp),  cos(pp)*cos(rr), 0,
		sin(pp), sin(rr)*cos(pp), 0, cos(pp)*cos(rr))

  P <- matrix(data = p_data, nrow = 4, ncol = 4)
  
  # Make resulting transformation matrix
  R <- H %*% P %*% trans_matrix1

  return(R)
}

xyz_2_enu <- function(trans_matrix1, res_trans, xyz){
  
  beam <- xyz_2_beam(trans_matrix1, xyz)
  
  # Given beam velocities, ENU coordinates are calculated as
  enu <- beam_2_enu(res_trans, beam) 
  return(enu)
}

enu_2_xyz <- function(trans_matrix1, res_trans, enu){
  
  xyz <- trans_matrix1 %*% solve(res_trans) %*% enu
   
  return(xyz)
}

beam_2_enu <- function(res_trans, beam){
  
  enu <- res_trans %*% beam
  
  return(enu)
}

beam_2_xyz <- function(trans_matrix1, beam){
  xyz <- trans_matrix1 %*% beam
}

xyz_2_beam <- function(trans_matrix1, xyz){
  
  beam <- solve(trans_matrix1) %*% xyz 
  return(beam)
}

enu_2_beam <- function(res_trans, enu){
  
  beam <- solve(res_trans) %*% enu 
  return(beam) 
}


########################################################################################
# Rotation of coordinates into mean flow direction  (mean(v)=0)
########################################################################################

# reference: Foken (2005)
# See Chapter 4 eqs 4.5, 4.6, 4.7 and 4.8 for this method.

project <- function(chunk.adv){

xyz <- chunk.adv[,1:3]

u=xyz[,1]
v=xyz[,2]
w=xyz[,3]

    A = atan2(mean(v),mean(u))
    
    urot = u*cos(A)+v*sin(A)
    vrot = -u*sin(A)+v*cos(A)
   
    umean = mean(urot)
    wmean = mean(w) 
    theta= atan ( wmean / umean ) 

    wcor = w * cos(theta) - urot * sin(theta) 
    ucor = w * sin(theta) + urot * cos(theta) 

result <- cbind(ucor,vrot,wcor)

return(result)
} 


#################################################
# along-stream mean velocity
#################################################

# code from GDopp package ( https://github.com/USGS-R/GDopp )

velocity_calc <- function(chunk.adv){
  
  veloc.cube <- chunk.adv[,1]^2 + chunk.adv[,2]^2 + chunk.adv[,3]^2
 
  veloc <- veloc.cube^(1/3)
  # delta time, will be divided off, so is arbitrary unless it is irregular
  del.t = 1 
  dist.traveled <- veloc*del.t
  adv.veloc <- sum(dist.traveled)/length(dist.traveled)
  
  return(adv.veloc) 
}


########################################################################################
# Despiking following Goring and Nikora (2002), and Mori et al. 2007
# R code translated by Marcus Klaus from Matlab code by Mori (2020)
########################################################################################

despikePhasespace3d <- function(fi, i_plot=TRUE, i_opt=2){
#======================================================================
#
# Version 1.2
#
# This subroutine excludes spike noise from Acoustic Doppler 
# Velocimetry (ADV) data using phase-space method, using 
# modified Goring and Nikora (2002) method by Nobuhito Mori (2005).
# Further modified by Joseph Ulanowski to remove offset in output (2014).
# 
#======================================================================
#
# Input
#   fi     : input data with dimension (n,1)
#   i_plot : =TRUE plot results (optional)
#   i_opt : = 0 or not specified  ; return spike noise as NaN
#           = 1            ; remove spike noise and variable becomes shorter than input length
#           = 2            ; interpolate NaN using cubic polynomial
#
# Output
#   fo     : output (filtered) data
#   ip     : excluded array element number in fi
#
# Example: 
#   [fo, ip] = func_despike_phasespace3d( fi, 9 );
#     or
#   [fo, ip] = func_despike_phasespace3d( fi, 9, 2 );
#
#
#======================================================================
# Terms:
#
#       Distributed under the terms of the terms of the BSD License
#
# Copyright:
#
#       Nobuhito Mori
#           Disaster Prevention Research Institue
#           Kyoto University
#           mori@oceanwave.jp
#
#========================================================================
#
# Update:
#       1.2     2014/03/18 Offset removed for non-zero mean data [J.U.]
#       1.11    2009/06/09 Interpolation has been added.
#       1.01    2009/06/09 Minor bug fixed
#       1.00    2005/01/12 Nobuhito Mori
#
#========================================================================

#
# --- initial setup
#
# number of maximum iternation
n_iter = 20
n_out  = 999
n      = length(fi)
f_mean = 0     # do not calculate f_mean here, as it will be affected by spikes (was: f_mean = nanmean(fi);)
f      = fi    # this offset subtraction is unnecessary now (was: f = fi - f_mean;)
lambda = sqrt(2*log(n))


#
# --- loop
#
n_loop = 1
while((n_out != 0) & (n_loop <= n_iter)){
#
# --- main
#
# step 0
f_mean=f_mean+mean(f,na.rm=TRUE) # accumulate offset value at each step [J.U.]
f = f - mean(f,na.rm=TRUE)
#nanstd[f]

# step 1: first and second derivatives
f_t  = gradient(f)
f_tt = gradient(f_t)

# step 2: estimate angle between f and f_tt axis
if(n_loop==1){
  theta = atan2( sum(f*f_tt), sum(f^2) ) 
}

# step 3: checking outlier in the 3D phase space
xp = func_excludeoutlier_ellipsoid3d(f,f_t,f_tt,theta)[[1]]
yp = func_excludeoutlier_ellipsoid3d(f,f_t,f_tt,theta)[[2]]
zp = func_excludeoutlier_ellipsoid3d(f,f_t,f_tt,theta)[[3]]
ip = func_excludeoutlier_ellipsoid3d(f,f_t,f_tt,theta)[[4]]
coefs = func_excludeoutlier_ellipsoid3d(f,f_t,f_tt,theta)[[5]]
#
# --- excluding data
#
n_nan_1 = length(which(is.na(f)==TRUE))
f[ip]  = NA
n_nan_2 = length(which(is.na(f)==TRUE))
n_out   = n_nan_2 - n_nan_1
#
# --- end of loop
#
n_loop = n_loop + 1
}


#
# --- post process
#
go = f + f_mean    # add offset back
ip = which(is.na(go))
if(n_loop < n_iter){
  print(
    paste(' Number of outlier   = ', num2str(length(which(is.na(f)==TRUE))), 
     ' : Number of iteration = ', num2str(n_loop-1))
  )
}else{
  print(
    paste(' Number of outlier   = ', num2str(length(which(is.na(f)==TRUE))), 
     ' : Number of iteration = ', num2str(n_loop-1), ' !!! exceed maximum value !!!' )
  )
}
#
# --- interpolation or shorten NaN data
#
if(abs(i_opt) >= 1){
	# remove NaN from data
    inan = which(!is.na(go))
    fo = go[inan]
    # interpolate NaN data
    if(abs(i_opt) == 2){
        


        x   = which(!is.na(go))
        y   = go[x]
        xi  = 1:length(fi)

        ## interp1 does not work if first or last value is NA. use NA for xi values outside this range 
        ## see https://stackoverflow.com/questions/47295879/using-interp1-in-r
        interp2 <- function(x, y, xi = x, ...) {
        yi <- rep(NA, length(xi));
        sel <- which(xi >= range(x)[1] & xi <= range(x)[2])
        yi[sel] <- interp1(x = x, y = y, xi = xi[sel], ...)
        return(yi)
        }

        fo = interp2(x, y, xi, 'cubic')
    }
}else{
    # output despiked value as NaN
    fo = go
}
#
# --- for check and  plot
#
if(i_plot == TRUE){ 
#theta/pi*180
F    = fi - f_mean
F_t  = gradient(F)
F_tt = gradient(F_t)
RF = rbind( 
   c(cos(theta),0,sin(theta)), 
   c(0,1,0), 
   c(-sin(theta),0,cos(theta)))

RB = rbind( 
   c(cos(theta),0,-sin(theta)), 
   c(0,1,0), 
   c(sin(theta),0,cos(theta)))


# making ellipsoid data
a = coefs[1]
b = coefs[2]
c = coefs[3]
ne  = 32
dt  = 2*pi/ne
dp  = pi/ne
t   = seq(from=0,to=2*pi,by=dt)
p   = seq(from=0,to=pi,by=dp)
n_t = length(t)
n_p = length(p)
# making ellipsoid
xe <- NULL
ye <- NULL
ze <- NULL
for(it in 1:n_t){
  for(is in 1:n_p){
    xe[n_p*(it-1)+is] = a*sin(p[is])*cos(t[is])
    ye[n_p*(it-1)+is] = b*sin(p[is])*sin(t[is])
    ze[n_p*(it-1)+is] = c*cos(p[is])
  }
}
xer = xe*RB[1,1] + ye*RB[1,2] + ze*RB[1,3]
yer = xe*RB[2,1] + ye*RB[2,2] + ze*RB[2,3]
zer = xe*RB[3,1] + ye*RB[3,2] + ze*RB[3,3]


# plot figures

## Fig 1
#cloud(f_tt~f*f_t,xlab="u",ylab="d u",zlab=("d^2 u"))
#  cloud(F_tt[ip]~F[ip]*F_t[ip],col=2,cex=1.2, add=TRUE)
#  cloud(zer~xer*yer,type="l", add=TRUE)

## Fig 2
plot(fi,type="l")
if(i_opt==2){
    points(fo,type="l",col=4)
}
points(ip,fi[ip],col=2)

} ## close (if plot yes)

## return filtered data and IDs of outliers
return(list(fo,ip))
}


func_excludeoutlier_ellipsoid3d <- function(xi,yi,zi,theta){
#======================================================================
#
# Version 1.01
#
# This program excludes the points outside of ellipsoid in two-
# dimensional domain
#
# Input
#   xi : input x data
#   yi : input y data
#   zi : input z data
#   theta  : angle between xi and zi
#
# Output
#   xp : excluded x data
#   yp : excluded y data
#   zp : excluded y data
#   ip : excluded array element number in xi and yi
#   coef : coefficients for ellipsoid
#
# Example: 
#   [xp,yp,zp,ip,coefs] = func_excludeoutlier_ellipsoid3d(f,f_t,f_tt,theta)
#
#
#======================================================================
# Terms:
#
#       Distributed under the terms of the terms of the BSD License
#
# Copyright: 
#
#       Nobuhito Mori, Kyoto University
#
#========================================================================
#
# Update:
#       1.01    2009/06/09 Nobuhito Mori
#       1.00    2005/01/12 Nobuhito Mori
#
#========================================================================
#
# --- initial setup
#
n = length(xi)
lambda = sqrt(2*log(n))
xp = NULL
yp = NULL
zp = NULL
ip = NULL
coefs = NULL
#
# --- rotate data
#
#theta = atan2( sum(xi.*zi), sum(xi.^2) )
if(theta == 0){
  X = xi
  Y = yi
  Z = zi
}else{
  R = rbind( 
   c(cos(theta),0,sin(theta)), 
   c(0,1,0), 
   c(-sin(theta),0,cos(theta)))

  X = xi*R[1,1] + yi*R[1,2] + zi*R[1,3]
  Y = xi*R[2,1] + yi*R[2,2] + zi*R[2,3]
  Z = xi*R[3,1] + yi*R[3,2] + zi*R[3,3]
}

#
# --- preprocess
#
a = lambda*sd(X,na.rm=T)
b = lambda*sd(Y,na.rm=T)
c = lambda*sd(Z,na.rm=T)
#
# --- main
#
m = 0
for (i in 1:n){
  x1 = X[i]
  y1 = Y[i]
  z1 = Z[i]
  # point on the ellipsoid
  x2 = a*b*c*x1/sqrt((a*c*y1)^2+b^2*(c^2*x1^2+a^2*z1^2))
  y2 = a*b*c*y1/sqrt((a*c*y1)^2+b^2*(c^2*x1^2+a^2*z1^2))
  zt = c^2* ( 1 - (x2/a)^2 - (y2/b)^2 )
  if( is.na(z1) ){
  z2 = 0
  } else {
  if (z1 < 0){
    z2 = -sqrt(zt)
  }else if ( z1 > 0){
    z2 = sqrt(zt)
  }else{
    z2 = 0
  }
  }
  # check outlier from ellipsoid
  dis = (x2^2+y2^2+z2^2) - (x1^2+y1^2+z1^2)
  if (is.na(dis)==FALSE){
  if (dis < 0){ 
    m = m + 1
    ip[m] = i
    xp[m] = xi[i]
    yp[m] = yi[i]
    zp[m] = zi[i]
  }
  }
}
coefs[1] = as.vector(a)
coefs[2] = as.vector(b)
coefs[3] = as.vector(c)

output <- list(xp,yp,zp,ip,coefs)
return(output)

}


################################################################################################
################################################################################################
################################################################################################
## Part 3: Data processing
################################################################################################
################################################################################################
################################################################################################

######################################
## Load data
######################################

## get ADV filenames
filenames <- substr(list.files(pattern="*.dat", full.names=TRUE),3,99) 
filenameshdr <- substr(list.files(pattern="*.hdr", full.names=TRUE),3,99)
 
## make empty output file
output <- NULL

######################################
## loop through all ADV files
######################################

for (i in 1:length(filenames)){

## extract metadata on ADV collection from the hdr file
DateTime <- gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=5)[5])) ## DateTime
VelocityRange <- gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=11)[11]))
VelocityRange  <- as.numeric(gsub("/", "",VelocityRange))
TransmitLength <- as.numeric(gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=12)[12])))
SamplingLength <- gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=13)[13]))
SamplingLength <- as.numeric(gsub("/", "",SamplingLength))
Distance <- gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=113)[113]))
Distance <- as.numeric(gsub("/", "",Distance))
Temperature <- gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=118)[118]))
Temperature <- as.numeric(gsub("/", "",Temperature))
VelocityScaling <- gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=17)[17]))
VelocityScaling <- as.numeric(gsub("/", "",VelocityScaling))
CoordinateSystem <- substr(readLines(paste0(folder.nm,"/",filenameshdr[i]),n=19)[19],39,42) ## coordinate system

## calculate water depth based on ADV distance-to-wall estimate
WaterDepth <- Distance*cos(position_data$pitch) + TransmitDepth

### calculate depth of sampling volume
### 5 cm distance between transmit transducer and sampling ("sensor") volume; add depth of transmit transducer
SensorDepth <- 0.05*cos(position_data$pitch) + TransmitDepth 

### extract ADV transformation matrix
Transformationmatrix1 <- as.numeric(unlist(strsplit((gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=82)[82])))," ")))
Transformationmatrix2 <- as.numeric(unlist(strsplit((gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=83)[83])))," ")))
Transformationmatrix3 <- as.numeric(unlist(strsplit((gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=84)[84])))," ")))
Transformationmatrix4 <- as.numeric(unlist(strsplit((gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=85)[85])))," ")))
trans_data <- c(Transformationmatrix1,Transformationmatrix2,Transformationmatrix3,Transformationmatrix4) 
trans_matrix <- matrix(data = trans_data, ncol = 4, byrow = TRUE)
# Scale the transformation matrix correctly to floating point numbers
trans_matrix <- trans_matrix/4096 

########################################
## load ADV velocity data file

chunk.adv <- load_adv(filenames[i],folder.nm)

#################################################################################################
### first despike and then transform coordinate system, following advice by Doroudian et al. (2010)

### despike data following Goring and Nikora (2002) with modifications by Mori et al. 2007

dataDespike <- chunk.adv[,3:6]
propDespiked <- matrix(0,1,4)
for(b in 1:4){
dataDespike[,b] <- despikePhasespace3d(chunk.adv[,b+2],plotting,2)[[1]]
propDespiked[,b] <- length(despikePhasespace3d(chunk.adv[,b+2],FALSE,2)[[2]])/length(dataDespike[,b])
local({
dev.set (2)
dev.print (device=tiff, file=paste0(filenames[i],substr(names(chunk.adv)[b+2],10,11),"_Timeseries.tiff"), width=par("din")*100, res=100);
})
dev.off ()
}

names(propDespiked) <- c("propDespiked.u","propDespiked.v","propDespiked.w1","propDespiked.w2")

## remove "NA" values at begin and end of time series, if present, to allow spectral plotting
dataDespike <- na.omit(dataDespike)

########################################
## transform coordinates (XYZ -> ENU)

dataDespike  <- coord_transform(trans_matrix, dataDespike, position_data)

#######################################################################################
## extract some useful velocity statistics and quality check measures from ADV data

## Signal-to-noise-ratio (SNR) + correlation 
SNR <- colMeans(chunk.adv[,11:14])
CORR <- colMeans(chunk.adv[,15:18])

## quality check: is Taylors hypothesis of frozen turbulence fulfilled?

## mean flow velocity in u direction
V <- velocity_calc(project(dataDespike))

FrozenTurb <- NULL
for(h in 1:4){
nrm.v. <- dataDespike[,h]-mean(dataDespike[,h]) ## fluctuating velocities
r.v. <- sqrt(sum(nrm.v.^2)/length(nrm.v.)) ## RMS of fluctuating velocities
FrozenTurb[h] <- (r.v./V)^3
} 
FrozenTurb <- matrix(FrozenTurb,ncol=4)
names(FrozenTurb) <- c("FrozenTurb.u","FrozenTurb.v","FrozenTurb.w1","FrozenTurb.w2")

## standard deviations of velocities
sdVel <-  apply(dataDespike[,1:4],2,FUN=sd)
names(sdVel) <- c("sdVel.u","sdVel.v","sdVel.w1","sdVel.w2")


################################ 
## calculate Turbulence indices

## Turbulent kinetic energy (m2/s2)
TKE <- TKEf(sdVel[1],sdVel[2],sdVel[3])
## Reynolds number
Re <- Ref(V=V,D=WaterDepth,kinvis=get_kin_viscosity(temperature = Temperature))
## Froude number
Fr <- Frf(V=V,D=WaterDepth)
## Energy dissipation rate (driven by drag) (Raymond et al. 2012)
eD <- Edf(S=Slope,V=V)
## Energy dissipation rate (driven by bed shear) (Raymond et al. 2012); assumes smooth bottom!
eS <- Esf(S=Slope,D=WaterDepth,W=0.4)


################################################
### calculate epsilon at sensor depth
################################################

Epsilon <- fit_epsilon(data=dataDespike,freq=200, WaterDepth=WaterDepth, SensorDepth=SensorDepth,SamplingLength=SamplingLength, diagnostic = plotting, filt=FALSE,n=12,fmax=10,SpAve=SpAve)
local({
dev.set (2)
dev.print (device=tiff, file=paste0(filenames[i],"_Spectra.tiff"), width=par("din")*100, res=100);
})
dev.off ()

############################
#### compile and save output
############################

output0 <- c(
filenames[i],
DateTime,
VelocityRange,
TransmitLength,
SamplingLength,
Distance,
Temperature,
VelocityScaling,
CoordinateSystem,
WaterDepth,
SensorDepth,
TKE,
Re,
Fr,
eD,
eS,
V,
SNR,
CORR,
FrozenTurb,
propDespiked,
sdVel,
unmatrix(Epsilon,byrow=TRUE))

names(output0)[c(1:17)] <- c("filename","DataTime","VelocityRange","TransmitLength","SamplingLength","Distance","Temperature",
"VelocityScaling","CoordinateSystem","WaterDepth","SensorDepth","TKE","Re","Fr","eD","eS","Vadv")

output <- rbind(output,output0)

# save updated model output dataframe
write.csv(output, file="output.csv")

} ## close i loop (through different experiments)


################################################################################################
################################################################################################
################################################################################################
## Part 4: Investigate and quality check output
################################################################################################
################################################################################################
################################################################################################

## read output file
output <- read.csv("output.csv")

###############################################################
### Post-processing quality check
###############################################################
### Note: different thresholds have been used in the literature
### I here use commonly used thresholds for SNR and correlations
###############################################################

## SNR should be > 15, correlation should be > 70 
SNRtest <- as.numeric(output[,which(colnames(output)=="correlation.X")]) > 70
Corrtest <- as.numeric(output[,which(colnames(output)=="signal.rat.X")]) > 15

## assumption of frozen turbulence should hold, i.e. (r.v./V)^3 < 1 (Kitaigorodskii et al. 1983)
Frozentest <- as.numeric(output[,which(colnames(output)=="FrozenTurb.u")]) < 1

## spectral fit should have positive R2
R2test <- as.numeric(output[,which(colnames(output)=="realr2.u")]) > 0

## extract means and 95% confidence intervals for epsilon estimates
epsUlwr <- as.numeric(output[,which(colnames(output)=="epsilonMLElwr.u")])
epsUupr <- as.numeric(output[,which(colnames(output)=="epsilonMLEupr.u")])
epsVlwr <- as.numeric(output[,which(colnames(output)=="epsilonMLElwr.v")])
epsVupr <- as.numeric(output[,which(colnames(output)=="epsilonMLEupr.v")])
epsW1lwr <- as.numeric(output[,which(colnames(output)=="epsilonMLElwr.w1")])
epsW1upr <- as.numeric(output[,which(colnames(output)=="epsilonMLEupr.w1")])
epsW2lwr <- as.numeric(output[,which(colnames(output)=="epsilonMLElwr.w2")])
epsW2upr <- as.numeric(output[,which(colnames(output)=="epsilonMLEupr.w2")])
epsUmean <- as.numeric(output[,which(colnames(output)=="epsilonMLE.u")])
epsVmean <- as.numeric(output[,which(colnames(output)=="epsilonMLE.v")])
epsW1mean <- as.numeric(output[,which(colnames(output)=="epsilonMLE.w1")])
epsW2mean <- as.numeric(output[,which(colnames(output)=="epsilonMLE.w2")])

## estimates for w1 and w2 should be similar, i.e. 95% confidence intervals should overlap (are "replicates")

Wtest <- epsW1lwr < epsW2upr &
epsW1upr > epsW2lwr 

## assumption of isotropic turbulence should hold, epsilon ratios for u,v,w1,w2 should be around 1
## check if 95% confidence intervals of estimates overlap for w1 compared to u and v
Isotest95 <- epsW1lwr < epsUupr &
epsW1upr > epsUlwr &
epsW1lwr < epsVupr &
epsW1upr > epsVlwr


##############################
## investigate results
############################## 

## investigate results of quality checks
data.frame(output[,which(colnames(output)=="filename")],SNRtest,Corrtest,Frozentest,R2test,Wtest,Isotest95)

## plot epsilon
plot(as.numeric(output[,which(colnames(output)=="epsilonMLE.u")]),ylim=c(0.00001,1),log="y")
points(as.numeric(output[,which(colnames(output)=="epsilonMLE.v")]),col=2)
points(as.numeric(output[,which(colnames(output)=="epsilonMLE.w1")]),col=3)
points(as.numeric(output[,which(colnames(output)=="epsilonMLE.w2")]),col=4)

## plot R2
plot(as.numeric(output[,which(colnames(output)=="realr2.u")]),ylim=c(-1,1))
points(as.numeric(output[,which(colnames(output)=="realr2.v")]),col=2)
points(as.numeric(output[,which(colnames(output)=="realr2.w1")]),col=3)
points(as.numeric(output[,which(colnames(output)=="realr2.w2")]),col=4)

######################################################################################################
# calculate geometric mean over "replicate" measurements for each experiment, following Baker and Gibson (1987)
######################################################################################################

## extract names of experimental runs from filenames (over which epsilon will be averaged)
exp <- substr(output[,which(colnames(output)=="filename")],1,14)

stats <- matrix(NA,length(unique(exp)),4)
for (j in 1:length(unique(exp))){
# select experiment and only observations that pass all quality checks
id <- which(exp==unique(exp)[j] & Isotest95==TRUE & SNRtest==TRUE & Corrtest==TRUE & Frozentest==TRUE & R2test==TRUE & Wtest)
## calculate mean and standard deviation of epsilon following Baker & Gibson (1987)
## epsilon (here "X") is log-distributed variable
X=as.numeric(output[id,which(colnames(output)=="epsilonMLE.w1")])
mu=mean(log(X),na.rm=T)
s2=var(log(X),na.rm=T)
## expected value of X
stats[j,1] <- exp(mu  + s2/2 )
# standard deviation of X
stats[j,2] <- sqrt(exp(2*mu + s2)*(exp(s2)-1))
# store number of observations that passed quality check
stats[j,3] <- length(id)
stats[j,4] <- length(which(exp==unique(exp)[j]))
}

## resulting means and sd of epsilon for each experiment
result <- data.frame(unique(exp),stats)
colnames(result) <- c("Experiment","EpsMean","EpsSD","nPass","n")

# save result
write.csv(result, file="Results/ExpMeansSDs.csv")






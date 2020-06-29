###############################################################################################################
### This R code calculates turbulent kinetic energy dissipation rate (epsilon) from 4D Nortek Vectrino ADV data
###############################################################################################################
### GENERAL INFORMATION
### code is written by Marcus Klaus, Swedish University of Agricultural Sciences (SLU), Umeå, Sweden
### contact: marcus.klaus@posteo.net
### last update: 2020-03-26
### 
### Please don't hesitate to contact me for any questions or if you find errors in the Code!
###
### Please note: There is no golden standard method for estimating epsilon in rivers
### results are sensitive to many parameters, and must always be quality checked!
###
### I also provide example data files ("ADV_ExampleFile") collected using a 4D Nortek Vectrino ADV 
### in a flume at Lunz Mesocosm facility, Lunz am See, Austria
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
### outgassing via Drifting and Anchored chambers:results of a ?ume experiment, for submission to Biogeosciences.
###  
### Code is based on methods in Zappa et al. 2003 and uses many improvements by Bluteau et al. 2011
### Some code parts are adapted from the GDopp package (https://github.com/USGS-R/GDopp)
###
### References:
###
### Bluteau, C., Jones, N., and Ivey, G.: Estimating turbulent kinetic energy dissipation using the inertial 
### subrange method in environmental ?ows, Limnology and Oceanography: Methods, 9, https://doi.org/10.4319/lom.2011.9.302, 2011.
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
###############################################################################################################
### OUTPUTS:  
# Variable name		Explanation
# filename			Filename
# DataTime			Date and time
# VelocityRange		Vectrino setting of Velocity range (m/s)
# TransmitLength		Vectrino setting of Transmit length (mm)
# SamplingVolume		Vectrino setting of Sampling volume (mm)
# Distance			Distance between sensor and "wall"
# Temperature		Water temperature (dC)
# VelocityScaling		Vectrino setting of Velocity scaling
# CoordinateSystem	Coordinate system of data aquisition
# WaterDepth		Water depth (m)
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
# decade:u			Ratio in power spectral density between upper and lower boundary of inertial subrange (u direction), based on unfiltered spectrum, 10 means the ratio is > 10; if ratio < 10, the number indicated shows the ratiom based on unfiltered spectrum
# decade:v			Ratio in power spectral density between upper and lower boundary of inertial subrange (v direction), based on unfiltered spectrum, 10 means the ratio is > 10; if ratio < 10, the number indicated shows the ratio
# decade:w1			Ratio in power spectral density between upper and lower boundary of inertial subrange (w1 direction), based on unfiltered spectrum, 10 means the ratio is > 10; if ratio < 10, the number indicated shows the ratio
# decade:w2			Ratio in power spectral density between upper and lower boundary of inertial subrange (w2 direction), based on unfiltered spectrum, 10 means the ratio is > 10; if ratio < 10, the number indicated shows the ratio
# k600.u			k600 (m/d), based on velocity time series in u direction; based on unfiltered spectrum
# k600.v			k600 (m/d), based on velocity time series in v direction; based on unfiltered spectrum
# k600.w1			k600 (m/d), based on velocity time series in w1 direction; based on unfiltered spectrum
# k600.w2			k600 (m/d), based on velocity time series in w2 direction; based on unfiltered spectrum

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


###############################################
## set directory to folder with ADV files
###############################################
 
setwd("C:/MARCUS/WORK/Applications/Aquacosm_TA/Data/ADV/")
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

###############################################
### load functions for hydrodynamic indices
###############################################

## turbulent kinetic energy dissipation [m2/s3] (Raymond et al. 2012)
#(driven by drag)
Edf <- function(S,V){
Ed=9.81*S*V
return(Ed)
}
#(driven by shear)
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


fit_epsilon <- function(dataset,freq=200, WaterDepth=WaterDepth, SensorDepth=SensorDepth,SamplingVolume=SamplingVolume,
diagnostic = FALSE,alpha=1.5*18/55*c(1,1.33,1.33,1.33),filt=TRUE,n=12,fmax=50){
  
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
  ## calculated as in Bluteau et al. 2011, using "project" function defined further down 
  v.mn <- velocity_calc(project(dataset))

  ##########################################  
  ## calculate characteristic length scales

  # conversion from frequency to wavenumber space (rad/m)
  wavenum <- 2*pi*wu$freq/v.mn 
  # length scale of sampling volume; divide by 1000 to convert mm to m
  vol.k = 2*pi*(1/SamplingVolume*1000) 
  # wave number corresponding to stream width
  width.k = 2*pi*(1/WaterWidth)
  # wave number corresponding to stream depth (calculated from ADV measurement)
  depth.k = 2*pi*(1/WaterDepth)
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
## fit spectrum for a series of candidate windows which where spectral density falls of by at least a decade  
## loop through candidate windows to find window with best model fit, i.e. lowest mean average deviation (MAD) (loosely following Bluteau et al. 2011)

## make output vectors
lower.k.s <- matrix(NA,1000,4)
upper.k.s <- matrix(NA,1000,4)
realr2.s <- matrix(NA,1000,4)
MAD.s <- matrix(NA,1000,4)
epsilonMLEmean.s <- matrix(NA,1000,4)
epsilonMLElwr.s <- matrix(NA,1000,4)
epsilonMLEupr.s <- matrix(NA,1000,4)
decade.s <- matrix(NA,1000,4)
epsilonMLE <- NULL
epsilonMLElwr <- NULL
epsilonMLEupr <- NULL
realr2 <- NULL
MAD <- NULL
lower.k <- NULL
upper.k <- NULL
decade <- NULL

## loop through all h ADV beam components
for (h in 1:4){

# set start values
s=2
# set upper window limit to wave number associated with ADV sampling volume
upper.k.s[1,h] <- vol.k
lower.k.s[1,h] <- depth.k

# loop through candidate windows
while (upper.k.s[s-1,h] <= vol.k){
lower.k.s[s,h] <- wavenum[ which(wavenum>= depth.k) ][s]

### find upper boundary, must be minimum a decade between upper and lower bound (Bluteau et al. 2011); 
### if no decade difference, reduce stepwise until upper boundary found, but do not accept drop to below half a decade
upper.k.s[s,h] <- wavenum[min(which(runmed(w.s[[h]],11)[which(lower.k.s[s,h] == wavenum)]/runmed(w.s[[h]],11) > 10  ))] 
decade.s[s,h]=10
dd=10
while(is.na(upper.k.s[s,h])){
if(dd<5){
break
}
upper.k.s[s,h] <- wavenum[min(which(runmed(w.s[[h]],11)[which(lower.k.s[s,h] == wavenum)]/runmed(w.s[[h]],11) > dd  ))] 
dd=dd-0.1
## store the ratio which was finally used
decade.s[s,h]=dd+0.1
}


## if there is no decade difference, stop loop
if(is.na(upper.k.s[s,h])==TRUE){
break
}
if(upper.k.s[s,h]> noise.k & filt==TRUE){
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
}
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
}
s

## select final (best) epsilon estimate and save estimate, goodness of fit measures and confidence interval
epsilonMLE[h] <- epsilonMLEmean.s[id,h]
epsilonMLElwr[h] <- epsilonMLElwr.s[id,h]
epsilonMLEupr[h] <- epsilonMLEupr.s[id,h]
realr2[h] <- realr2.s[id,h]
MAD[h] <- MAD.s[id,h]
lower.k[h] <- lower.k.s[id,h]
upper.k[h] <- upper.k.s[id,h]
decade[h] <- decade.s[id,h]
} ## close h loop

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
    #abline(v=lower.k,col=c(1:4))
    #text(lower.k*1.15,5*10^-3,"lower",srt=90)
    #abline(v=upper.k,col=c(1:4))
    #text(upper.k*1.15,5*10^-3,"upper",srt=90)
    abline(v=width.k,lty=2)
    text(width.k*1.15,5*10^-3,"width",srt=90)
    abline(v=depth.k,lty=2)
    text(depth.k*1.15,5*10^-3,"depth",srt=90)
    abline(v=depth.k,lty=2)
    #abline(v=noise.k,lty=3)


    ## add lines based on which epsilon is estimated
    for (p in 1:4){
    ip <- wavenum>lower.k[p] & wavenum<upper.k[p]
    lines(wavenum[ip],alpha[p]*(epsilonMLE[p]^(2/3))*wavenum[ip]^(-5/3),col=p+4,lwd=2)
    # add confidence intervals
    #lines(wavenum[ip],alpha[p]*(epsilonMLElwr[p]^(2/3))*wavenum[ip]^(-5/3),col=p+4,lwd=1)
    #lines(wavenum[ip],alpha[p]*(epsilonMLEupr[p]^(2/3))*wavenum[ip]^(-5/3),col=p+4,lwd=1)
    }

    text(5,3*10^-7,paste0("R2u=",round(realr2[1],2)))
    text(5,10^-7,paste0("R2v=",round(realr2[2],2)))
    text(5,3*10^-8,paste0("R2w1=",round(realr2[3],2)))
    text(5,10^-8,paste0("R2w2=",round(realr2[4],2)))

    text(5,3*10^-5,paste0("eu=",round(epsilonMLE[1],4)))
    text(5,10^-5,paste0("ev=",round(epsilonMLE[2],4)))
    text(5,3*10^-6,paste0("ew1=",round(epsilonMLE[3],4)))
    text(5,10^-6,paste0("ew2=",round(epsilonMLE[4],4)))


### kolmogorov length scale (threshold for noise)
k.vis <- get_kin_viscosity(Temperature)
nk <- (k.vis^3/epsilonMLE)^(1/4)
nk.k <- 1/(10*nk)
abline(v=nk.k,lty=3,col=c(1:4))

legend("topright",c("u","v","w1","w2"),col=c(1,2,3,4),lwd=1,bg="white")

}

## frequency associated with kolmogorov scale; here I use it as input for filtering the spectrum
f= nk.k*v.mn/(2*pi)

## generate output
out <- rbind(epsilonMLE,epsilonMLElwr,epsilonMLEupr,realr2,MAD,lower.k,upper.k,f,decade)
colnames(out) <- c("u","v","w1","w2")
return(out)
}


###########################################################
# Function to convert epsilon to k600
###########################################################
# based on Zappa et al. (2003)
# code from GDopp package ( https://github.com/USGS-R/GDopp )

epsilon_to_k <- function(epsilon,temperature=20,nu=0.2){
  
  if (length(epsilon) != length(temperature)){
    stop('input vectors for epsilon and temperature must have the same number of elements')
  }
  k.vis <- get_kin_viscosity(temperature)
  m4s4 <- 86400^4
  e.k <- epsilon*k.vis*m4s4 #now in m/day
  
  # convert to m/day
  k600 <- nu*(e.k)^0.25*600^(-0.5)
  return(k600)
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
# I received 4*4 tilt and heading matrices by personal communication with Nortek on 2020-01-17

# transform velocity data between beam coordinates (along ADV beam direction) and ENU coordinates, where
# E = East-West component, N = North-South component, U = Up-Down component

# transformation follows
# Lohrmann, Atle, Ramon Cabrera, and Nicholas C. Kraus: Acoustic-Doppler velocimeter (ADV) for laboratory use. 
# Fundamentals and advancements in hydraulic measurements and experimentation, pp. 351-365. ASCE, 1994.


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

# reference: Thomas Foken (2005): Micrometeorology. Srpinger.
# https://link.springer.com/book/10.1007%2F978-3-540-74666-9
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

##########################################################################
## loop through all ADV files
## Note: this script assumes the same boundary conditions for all files
## adjust, if needed!
##########################################################################

for (i in 1:length(filenames)){

## extract metadata on ADV collection from the hdr file
DateTime <- gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=5)[5])) ## DateTime
VelocityRange <- gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=11)[11]))
VelocityRange  <- as.numeric(gsub("/", "",VelocityRange))
TransmitLength <- as.numeric(gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=12)[12])))
SamplingVolume <- gsub("[[:alpha:]]", "", sub("\\D+","",readLines(paste0(folder.nm,"/",filenameshdr[i]),n=13)[13]))
SamplingVolume <- as.numeric(gsub("/", "",SamplingVolume))
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
trans_data <- c(Transformationmatrix1,Transformationmatrix2,Transformationmatrix3,Transformationmatrix4) # trans matrix for example package data
trans_matrix <- matrix(data = trans_data, ncol = 4, byrow = TRUE)
# Scale the transformation matrix correctly to floating point numbers
trans_matrix <- trans_matrix/4096 

########################################
## load ADV velocity data file

chunk.adv <- load_adv(filenames[i],folder.nm)
chunk.advOrig <- chunk.adv[,3:ncol(chunk.adv)]

########################################
## transform coordinates (XYZ -> ENU)

chunk.adv <- coord_transform(trans_matrix, chunk.adv, position_data)
data <- chunk.adv

#################################################################################################
### despike data (using despike function from "oce" package, following Goring and Nikora 2002
### note: adjust the n and k parameters of the filter to your data!

dataDespike <- data
propDespiked <- matrix(0,1,4)
for(b in 1:4){
dataDespike[,b] <- despike(data[,b], reference = c("median"), n = 1.5, k = 41, min = NA, max = NA, replace = c("NA"))
propDespiked[,b] <- length(which(is.na(dataDespike[,b])))/length(dataDespike[,b])
dataDespike[,b] <- despike(data[,b], reference = c("median"), n = 1.5, k = 41, min = NA, max = NA, replace = c("reference"))
}

names(propDespiked) <- c("propDespiked.u","propDespiked.v","propDespiked.w1","propDespiked.w2")

########################################
## plot raw and despiked time series
########################################
## note: adjust plot limits, if needed
plot(ts(data[,1], frequency=200),type="l",col=5,ylab="Velocity [m/s]",xlab="Time [s]",ylim=c(-0.3,0.3))
points(ts(data[,2], frequency=200),type="l",col=6)
points(ts(data[,3], frequency=200),type="l",col=7)
points(ts(data[,4], frequency=200),type="l",col=8)
points(ts(dataDespike[,1], frequency=200),type="l",col=1)
points(ts(dataDespike[,2], frequency=200),type="l",col=2)
points(ts(dataDespike[,3], frequency=200),type="l",col=3)
points(ts(dataDespike[,4], frequency=200),type="l",col=4)

## save plots
legend("topleft",col=c(1,2,3,4),lwd=1,c("North","East","Up1","Up2"),ncol=4)
local({
dev.set (2)
dev.print (device=tiff, file=paste0(filenames[i],"_Timeseries.tiff"), width=par("din")*100, res=100);
})
dev.off ()

#######################################################################################
## extract some useful velocity statistics and quality measures from ADV data

## Signal-to-noise-ratio (SNR) + correlation 
SNR <- colMeans(chunk.advOrig[,9:12])
CORR <- colMeans(chunk.advOrig[,13:16])


## quality check: is Taylors hypothesis of frozen turbulence fulfilled?

FrozenTurb <- NULL
for(h in 1:4){
nrm.v. <- dataDespike[,h]-mean(dataDespike[,h]) ## fluctuating velocities
r.v. <- sqrt(sum(nrm.v.^2)/length(nrm.v.)) ## RMS of fluctuating velocities
FrozenTurb[h] <- (r.v./velocity_calc(project(dataDespike)))^3
} 
FrozenTurb <- matrix(FrozenTurb,ncol=4)
names(FrozenTurb) <- c("FrozenTurb.u","FrozenTurb.v","FrozenTurb.w1","FrozenTurb.w2")

## standard deviations of velocities
sdVel <-  apply(dataDespike[,1:4],2,FUN=sd)
names(sdVel) <- c("sdVel.u","sdVel.v","sdVel.w1","sdVel.w2")

## mean flow velocity in u direction
V <- abs(mean(project(dataDespike)[,1]))


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
## Energy dissipation rate (driven by bed shear) (Raymond et al. 2012)
eS <- Esf(S=Slope,D=WaterDepth,W=0.4)

################################################
### calculate epsilon
################################################
## should spectra be plotted?
dia = TRUE

## calculate for unfiltered spectrum
Epsilon <- fit_epsilon(data=dataDespike,freq=200, WaterDepth=WaterDepth, SensorDepth=SensorDepth,SamplingVolume=SamplingVolume, diagnostic = dia, filt=FALSE,n=12,fmax=10)
local({
dev.set (2)
dev.print (device=tiff, file=paste0("Results/",filenames[i],"_",substr(names(chunk.adv)[b],10,12),"_Spectra.tiff"), width=par("din")*100, res=100);
})
dev.off ()

## calculate for filtered spectrum (not used here to provide any output
#EpsilonFilt <- fit_epsilon(data=dataDespike,freq=200, WaterDepth=WaterDepth, SensorDepth=SensorDepth,SamplingVolume=SamplingVolume, diagnostic = dia, filt=TRUE,n=1,fmax=Epsilon[9,])
#local({
#dev.set (2)
#dev.print (device=tiff, file=paste0("Results/",filenames[i],"_",substr(names(chunk.adv)[b],10,12),"_Spectra_Filt.tiff"), width=par("din")*100, res=100);
#})
#dev.off ()
#rownames(EpsilonFilt) <- paste0(rownames(EpsilonFilt),"_filt")

# calculate k600 (m/d); "nu" scaling coefficient set to 0.16 following Moog and Jirka 1999 (for low-gradient gravel flume)
# note: "nu" coefficient is highly system-specific!

k600 <- matrix(0,1,4)
for(h in 1:4){
k600[,h] <- epsilon_to_k(Epsilon[2,h],temperature=Temperature,nu=0.16) 
}
names(k600) <- c("k600.u","k600.v","k600.w1","k600.w2")


############################
#### compile and save output
############################

output0 <- c(
filenames[i],
DateTime,
VelocityRange,
TransmitLength,
SamplingVolume,
Distance,
Temperature,
VelocityScaling,
CoordinateSystem,
WaterDepth,
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
unmatrix(Epsilon,byrow=TRUE),
k600)

names(output0)[c(1:16)] <- c("filename","DataTime","VelocityRange","TransmitLength","SamplingVolume","Distance","Temperature",
"VelocityScaling","CoordinateSystem","WaterDepth","TKE","Re","Fr","eD","eS","Vadv")

output <- rbind(output,output0)

# save updated model output dataframe
write.csv(output, file="Results/output.csv")

} ## close i loop (through different experiments)

## read output file
output <- read.table("output.txt", header=TRUE,sep = "\t")


###############################################################
### Post-processing quality check
###############################################################
### Note: different thresholds have been used in the literature
### I here use commonly used thresholds for SNR and correlations
###############################################################

## SNR should be > 15, correlation should be > 70 
as.numeric(output[,which(colnames(output)=="correlation.X")]) > 70
as.numeric(output[,which(colnames(output)=="signal.rat.X")]) > 15

## assumption of frozen turbulence should hold, i.e. (r.v./V)^3 < 1 (Kitaigorodskii et al. 1983)
as.numeric(output[,which(colnames(output)=="FrozenTurb.u")]) < 1

## spectral fit should have positive R2
as.numeric(output[,which(colnames(output)=="realr2:u")]) > 0


## check similarity in estimates for different directions
epsUlwr <- as.numeric(output[,which(colnames(output)=="epsilonMLElwr:u")])
epsUupr <- as.numeric(output[,which(colnames(output)=="epsilonMLEupr:u")])
epsVlwr <- as.numeric(output[,which(colnames(output)=="epsilonMLElwr:v")])
epsVupr <- as.numeric(output[,which(colnames(output)=="epsilonMLEupr:v")])
epsW1lwr <- as.numeric(output[,which(colnames(output)=="epsilonMLElwr:w1")])
epsW1upr <- as.numeric(output[,which(colnames(output)=="epsilonMLEupr:w1")])
epsW2lwr <- as.numeric(output[,which(colnames(output)=="epsilonMLElwr:w2")])
epsW2upr <- as.numeric(output[,which(colnames(output)=="epsilonMLEupr:w2")])

## estimates for w1 and w2 should be similar, i.e. 95% confidence intervals should overlap (are "replicates")

epsW1lwr < epsW2upr 
epsW1upr > epsW2lwr 


## assumption of isotropic turbulence should hold, epsilon ratios for u,v,w1,w2 should be around 1
epsUW1 <- as.numeric(output[,which(colnames(output)=="epsilonMLE_filt:u")])/
as.numeric(output[,which(colnames(output)=="epsilonMLE_filt:w1")]) 
epsUW2 <- as.numeric(output[,which(colnames(output)=="epsilonMLE_filt:u")])/
as.numeric(output[,which(colnames(output)=="epsilonMLE_filt:w2")]) 
epsVW1 <- as.numeric(output[,which(colnames(output)=="epsilonMLE_filt:v")])/
as.numeric(output[,which(colnames(output)=="epsilonMLE_filt:w1")]) 
epsVW2 <- as.numeric(output[,which(colnames(output)=="epsilonMLE_filt:v")])/
as.numeric(output[,which(colnames(output)=="epsilonMLE_filt:w2")]) 
epsW1W2 <- as.numeric(output[,which(colnames(output)=="epsilonMLE_filt:w1")])/
as.numeric(output[,which(colnames(output)=="epsilonMLE_filt:w2")]) 

## check if 95% confidence intervals of estimates overlap for w1 compared to u and v
epsW1lwr < epsUupr &
epsW1upr > epsUlwr &
epsW1lwr < epsVupr &
epsW1upr > epsVlwr

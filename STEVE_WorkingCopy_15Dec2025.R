# Rothamsted Park Grass parameters Plot 92a (limed to pH 7)
# Data for As

# Global variables --------------------------------------------------------

## General parameters
iTime <- 33 #52 # Simulation time (yrs)
iSimulationDepth <- 70 # Simulation depth (cm)
iIncorporationDepth <- 1  # The depth (cm) to which CTE is incorporated into the soil. Must be >0. If not incorporated, TE is still added to top cm

## Climate parameters
iRainfall <- Rain1907_1939$TotalRF # Rainfall and irrigation (mm/yr). Infiltrates and is taken up by plants down to root depth
iET <- 448 # Evapotranspiration (mm/yr). Taken up evenly by plants over root depth - despite uneven root distribution

### Yearly rainfall 
# This function allows iRainfall to be either a single annual average value or a vector of yearly values. 
# If yearly data is available, set iRainfall to a vector of yearly values that is equal in length to the value of iTime.
if (length(iRainfall) == 1) {
  RainfallSeries <- rep(iRainfall, iTime)
} else if (length(iRainfall) == iTime) {
  RainfallSeries <- iRainfall
} else {
  stop("Length of iRainfall must be 1 or equal to iTime")
}

## Soil parameters 
iSoilDensity <- 1.1 # Bulk density (t/m^3)
iSliceWaterMass <- 50 #T he mass of water (t/ha) that a 1cm slice of soil can hold. Essentially the porosity (%)
iSoilLoss <- 0.72 # The annual soil loss (t/ha)
iInitialSoilCon <- 13.02 # 23.82 # Initial concentration of TE in soil (mg/kg). Must be > 0
iKdSoil <- 440 # Soil adsorption coefficient

### Kd parameters
iUseLangmuir <- FALSE # switch: FALSE = constant Kd, TRUE = Langmuir. If FALSE, no need to specify iLang_Qmax and iLang_K
iLang_Qmax <- 23 # Langmuir maximum sorption capacity (same units as soil conc × mass)
iLang_K <- 0.01 # Langmuir affinity constant (same units⁻¹)

## Plant parameters
iRootDepth <- 36 # The root depth (cm). Must be <= simulation depth
iDryBiomassProduced <- 4.684 # 5.374 # Yield of biomass produced (t/ha)
iPercentRemoved <- 100 # Percentage of biomass removed
iK <- 0 # Decay constant (must be between 0-0.5, applicable to TEs that have plant-regulated uptake, particularly Cu. If not paplicable, set to 0)
iInitialCropCon <- 0.439 # 0.100 # Initial crop conc. (mg/kg d.w.) - must be >0 otherwise RAF will be 0

## Contaminant parameters
iConApplied <- 373 # 4.0 #annual amount of TE (g) that is added per hectare
iAtmosDeposition <- 2.5 # Atmospheric deposition (g/ha/yr)
iConAdded <- iAtmosDeposition + iConApplied

# Output variables --------------------------------------------------------

SoilCon <<- c(rep(0,iTime))
PlantData <<- c(rep(0,iTime))
Leached <<- c(rep(0,iTime))
LeachedConc <<- c(rep(0,iTime))
FinalProfile <<- c(rep(0,iSimulationDepth))
YearlyProfile <<- c(rep(0,iSimulationDepth))
top23 <<- c(rep(0,iTime))

# Mass balance variables --------------------------------------------------

MassAdded <- 0
MassInPlants <- 0
MassLeached <- 0
MassRunoff <- 0
MassInSoilStart <- 0
MassInSoilFinish <- 0


# Setting up slices -------------------------------------------------------

# soil is treated as every cm down to the simulation depth being its own 'slice' or layer
SoilMass <- iSoilDensity*100 # S[d]
ConMass <- SoilMass*iInitialSoilCon # M[d]

for (i in 1:iSimulationDepth)
{
  # setting up slices with the Contam mass (grams)
  SoilMass[i] <- iSoilDensity * 100 # soil mass for every slice
  TotalSoilMass <- sum(SoilMass) # The total mass of soil in the system - used to calculate average concentration
  ConMass[i] <- SoilMass[i] * iInitialSoilCon # contaminant mass in every slice
  MassInSoilStart <- sum(ConMass) # initial mass in soil - mass balance variable
}
# above code gives every slice appropriate g of contaminant based on initial soil conc


# Kd function ———---------------------------------------------------------

getKd <- function(C_soil) {
  if (!iUseLangmuir) {
    iKdSoil
  } else {
    # Langmuir: q = (Qmax * K * C) / (1 + K * C) → Kd = q/C
    (iLang_Qmax * iLang_K) / (1 + iLang_K * C_soil)
  }
}

# Solution Mass Function --------------------------------------------------

# Used in drainage and plant uptake
#SolutionMass <- function(i)
#{
#  ((ConMass[i] / SoilMass[i]) / iKdSoil) * iSliceWaterMass # (g)
#}

#RAF1 <- iInitialCropCon / (iInitialSoilCon / iKdSoil) # root adsorption factor -   initialized here and given a value outside


# New function using Langmuir Kd
SolutionMass <- function(i) {
  C_soil  <- ConMass[i] / SoilMass[i]
  Kd_curr <- getKd(C_soil)
  (C_soil / Kd_curr) * iSliceWaterMass
}

# Existing ISC calc
InitSolnConc <- SolutionMass(1) / iSliceWaterMass

# New functions for Initial RAF
Kd_init <- getKd(iInitialSoilCon)
RAF1 <- iInitialCropCon / (iInitialSoilCon / Kd_init)


# Incorporate -------------------------------------------------------------

Incorporate <- function(ContaminantMass,IncorporationDepth)
{
  SliceConAdded <- rep(c(ContaminantMass / IncorporationDepth, 0),
                       c(IncorporationDepth, iSimulationDepth-IncorporationDepth)) 
  ConMass <<- ConMass + SliceConAdded
}

# Plant Uptake ------------------------------------------------------------

PlantUptake <- function(HoldX)
{
  BiomassRemoved <- iPercentRemoved / 100 # the fraction of biomass that is removed from the site
  PlantConMass <- 0
  RemovedFromSlice <- 0
  for (i in 1: iRootDepth)
  {
    #SolnConc <- SolutionMass(i)/iSliceWaterMass
    
    # New function 2 May 2025 for Langmuir Kd
    C_soil   <- ConMass[i] / SoilMass[i]
    Kd_curr  <- getKd(C_soil)
    SolnConc <- (C_soil / Kd_curr)  # g contaminant per t water
    
    RAF <- (RAF1 * InitSolnConc) / (InitSolnConc + iK * (SolnConc - InitSolnConc)) # root adsorption factor - initialized here and given a value outside
    RemovedFromSlice <- (SolnConc * RAF * iDryBiomassProduced) / iRootDepth
    PlantConMass <- PlantConMass + RemovedFromSlice
    ConMass[i] <<- ConMass[i] - RemovedFromSlice
  }
  
  PlantData[HoldX] <<- (PlantConMass / iDryBiomassProduced)
  MassInPlants <<- MassInPlants + PlantConMass * BiomassRemoved # The fraction of biomass that is removed from the site
  Incorporate(PlantConMass * (1 - BiomassRemoved), iIncorporationDepth) # Reincorporation of the unused plant biomass - Updates ConMass again
}

# Drainage & Leaching ----------------------------------------------------

# New function to test:

Drainage <- function(HoldX)
{
  Leached <- 0
  MobileCon <- 0
  RainfallMass <- 10 * pmax(RainfallSeries[HoldX] - iET, 0) # Mass of water through soil profile (t/ha)
  
  Loops <- trunc(RainfallMass / iSliceWaterMass) # Number of passes rainfall makes through each slice
  
  for (steps in 1:Loops)
  {
    for (i in 1:iSimulationDepth)
    {
      SLNMS <- SolutionMass(i) # Solution mass for slice
      ConMass[i] <<- ConMass[i] - SLNMS
      
      if (i < iSimulationDepth)
      {
        ConMass[i + 1] <<- ConMass[i + 1] + SLNMS
      }
      else
      {
        Leached <- Leached + SLNMS
      }
    }
  }
  Leached[HoldX] <<- Leached
  MassLeached <<- MassLeached + Leached
  
  # Calculate leachate concentration - THERE IS AN ERROR HERE, it is likely now fixed 
  LeachedConc[HoldX] <<- Leached / RainfallMass # - (iSimulationDepth * iSliceWaterMass))) 
}


## Old function:

#Drainage <- function(HoldX)
#{
#  Leached <- 0
#  MobileCon <-0
#  RainfallMass <- 10*pmax(iRainfall-iET,0) #Mass of water going through the soil profile (t/ha)
  
#  Loops <- trunc(RainfallMass / iSliceWaterMass) # How many times rainfall will go thru a slice
  
#  for (steps in 1: Loops)
#  {
    
#    for (i in 1: iSimulationDepth)
#    {
#      SLNMS <- SolutionMass(i) #Need to assign this once so that solution mass is not recalculated half way through the tipping bucket
#      ConMass[i] <<- ConMass[i] - SLNMS
      
#      if (i < iSimulationDepth)
#      {
#        ConMass[i + 1] <<- ConMass[i + 1] + SLNMS
#      }
#      else
#      {
#        Leached <- Leached + SLNMS
#      }
#      
#    }
#  }
#  Leached[HoldX] <<- Leached
#  MassLeached <<- MassLeached + Leached
  
  # New line 10/6/24 to get leached concentration
#  LeachedConc[HoldX] <<- (Leached / (RainfallMass - (iSimulationDepth * iSliceWaterMass))) 
#}

# Soil Loss ---------------------------------------------------------------

SoilLoss <- function()
{
  SoilLossFract <- pmin(iSoilLoss / SoilMass[1],1)
  
  SoilConLoss <- ConMass[1] * SoilLossFract # This just removes contaminant from the 1st slice, depending on the amount of soil loss (runoff)
  ConMass[1] <<- ConMass[1] - SoilConLoss
  MassRunoff <<- MassRunoff + SoilConLoss
}

# Main Loop ---------------------------------------------------------------

#MainLoop <- function(iRootDepth, iRainfall, iET, iSoilDensity, iKdSoil, 
#                     iConApplied, iSoilLoss, iInitialCropCon, iInitialSoilCon, 
#                     iDryBiomassProduced, iPercentRemoved, iSliceWaterMass, iAtmosDeposition) 
MainLoop <- function() 
{
  for (year in 1: iTime)
  {
    #Incorporation
    MassAdded <<- MassAdded + iConAdded
    Incorporate(iConAdded, iIncorporationDepth) # adding the contaminant annually
    
    SoilLoss()
    Drainage(year)
    PlantUptake(year)
    
    #annual mass sum
    Mass <- sum(ConMass)
    for (i in 1: iSimulationDepth)
    {
      SoilCon[year] <<- (Mass / TotalSoilMass) # calculating the total Con mass (grams) in soil
    }
    #yearly profile
    for (i in 1: iSimulationDepth)
    {
      YearlyProfile[i] <<- (ConMass[i] / SoilMass[i])
    }
    top23[year] <<- mean(YearlyProfile[1:23])
  }
  MassInSoilFinish <<- sum(ConMass)
  for (i in 1: iSimulationDepth)
  {
    FinalProfile[i] <<- (ConMass[i] / SoilMass[i])
  }
  
  checksum <<- (MassInSoilFinish + MassLeached + MassInPlants + MassRunoff) - (MassInSoilStart + MassAdded) # checking the mass balance adds up
  top23Final <<- mean(FinalProfile[1:23])
  
}

MainLoop()


# Outputs -----------------------------------------------------------------

## Mass balance variables
print(paste ("initial mass in soil ",  MassInSoilStart))
print(paste ("mass added ", MassAdded))
print(paste ("mass in plants ",MassInPlants))
print(paste ("mass leached ",  MassLeached))
print(paste ("mass in runoff ", MassRunoff))
print(paste ("final mass in soil ", MassInSoilFinish))
print(paste ("checksum ",checksum))
print(paste ("topsoil concentration", top23Final))

# Mass balance check
if (abs(checksum) <0.001)
{
  print("OKAY")
}else
{
  print("error")
}

# basic plots
par(mfrow=c(2,2))
plot(FinalProfile, 1:iSimulationDepth, ylim = rev(range(1:iSimulationDepth)))
plot(1:iTime, PlantData)
plot(1:iTime, Leached)
plot(1:iTime, SoilCon)

# Data export for output plots --------------------------------------------
library(tidyverse)

DataExport <- data.frame(1:iSimulationDepth, FinalProfile)
MBExport <- data.frame(MassInSoilStart, MassAdded, MassRunoff, MassLeached, MassInPlants, MassInSoilFinish, checksum, top23Final)
ValidationExport <- data.frame(1:iTime, top23)
# ^ Need to adjust the above to be the concentration to the appropriate depth to compare to the validation data 

#dir.create(paste0("Outputs2025"), showWarnings = FALSE) #stops warnings if folder already exists
#write.csv(DataExport, row.names = F, file.path(paste0("Outputs2025"), "profileconc1947.csv"))
#write.csv(MBExport, row.names = F, file.path(paste0("Outputs2025"), "mass_bal1947.csv"))
#write.csv(ValidationExport, row.names = F, file.path(paste0("Outputs2025"), "Validation_ConcTime1947.csv"))


SoilConcs1 <- data.frame(1:iTime, top23)
#SoilConcs2 <- data.frame(1:iTime, top23)
SoilConcs <- bind_rows(SoilConcs1, SoilConcs2)
write.csv(SoilConcs, row.names = F, file.path(paste0("Outputs2025"), 'SoilConc.csv'))

PlantConcs1 <- data.frame(1:iTime, PlantData)
#PlantConcs2 <- data.frame(1:iTime, PlantData)
PlantConcs <- bind_rows(PlantConcs1, PlantConcs2)
#write.csv(PlantConcs, row.names = F, file.path(paste0("Outputs2025"), 'PlantConc.csv'))

LeachedConcs1 <- data.frame(1:iTime, LeachedConc)
#LeachedConcs2 <- data.frame(1:iTime, LeachedConc)
LeachedConcs <- bind_rows(LeachedConcs1, LeachedConcs2)
#write.csv(LeachedConcs, row.names = F, file.path(paste0("Outputs2025"), 'LeachedConc.csv'))


# Validation concentrations - no longer needed -----------------------------------------------


#i10 <- data.frame(iTime, top23)
#i28 <- data.frame(iTime, top23)
#i30 <- data.frame(iTime, top23)
#i47 <- data.frame(iTime, top23)
#i63 <- data.frame(iTime, top23)

#PlantConcs1 <- data.frame(1:iTime, PlantData)

year <- 1971
i12.2 <- data.frame(iTime, top23)
#i17.2 <- data.frame(iTime, top23)
#si25.2 <- data.frame(iTime, top23)
#i32.2 <- data.frame(iTime, top23, year)
#i136.2 <- data.frame(iTime, top23, year)
#i139.2 <- data.frame(iTime, top23, year)
#i43.2 <- data.frame(iTime, top23, year)
#i46.2 <- data.frame(iTime, top23, year)
#i49.2 <- data.frame(iTime, top23, year)
#i52.2 <- data.frame(iTime, top23, year)

#PlantConcs2 <- data.frame(1:iTime, PlantData)

library(tidyverse)
PlantConcs <- bind_rows(PlantConcs1, PlantConcs2)
#write.csv(PlantConcs, 'PlantConc.csv')

SoilConcs <- bind_rows(i10, i28, i30, i47, i63, i12.2, i17.2, si25.2, i32.2, i136.2, i139.2, i43.2, i46.2, i49.2, i52.2)
#write.csv(SoilConcs, "SoilConcs.csv")


SoilConcs2 <- data.frame(1:iTime, top23)



# Version info ------------------------------------------------------------

# Date      | Name         | Changes
---------------------------------------------------------------------------
# 6/6/2024  | H Thompson  | Changed main loop to give read out of topsoil concentration for every year of the simulation
# 10/6/24   | H Thompson  | Wrote in line to calculate concentration of TE in leachate at end of simulation 
# Oct 2025 Functions for Langmuir Kd and yearly weather data integrated

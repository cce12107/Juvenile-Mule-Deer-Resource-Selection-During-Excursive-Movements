# Early-life movement ecology of a cervid: Implications for Chronic Wasting Disease Management Analysis Code
# Calvin C. Ellis
# 12/17/2024

set.seed(4432)

##################################
# Identifying Excursive Movements
##################################
library(sf)
library(dplyr)
library(data.table)

for (i in 1:length(id)) { 
  ind_gps = gps_sf[gps_sf$id == id[i],]
  whichID <- id[i]
  
  #read in shapefiles
  natal = st_read(file.path(natalpath, paste(whichID,"_Natal.shp", sep="" )))
  
  #create intersection of points and natal range
  out <- st_intersection(ind_gps, natal)
  
  #create buffer dist based on steplength inside natal range
  bufferdist = mean(out$sl_)
  
  buffered_natal <- st_buffer(natal, dist=bufferdist)
  
  #find outside buffer points using st_intersect
  outsidepoints = st_intersects(ind_gps, buffered_natal)
  
  #create true/false list for if outside buffer or not
  outsidepoints = sapply(st_intersects(ind_gps, buffered_natal),function(x){length(x)==0})
  
  #create a column in gps dataframe based on true/false list
  ind_gps$outsidebuffer = outsidepoints
  
  #create column counting consecutive true/false outside buffers to determine which true's are greater than or equal to 12 consecutive points (24 hours)
  ind_gps %>% group_by(gr = cumsum(outsidebuffer != lag(outsidebuffer, default = first(outsidebuffer)))) %>% mutate(count = n())
  
  setDT(ind_gps)[, count:= .N, rleid(outsidebuffer)]
  
  #create a column for excursion or not (y or n) based on outside buffer and consecutive gps points
  ind_gps$excursion = ifelse(
    (
      (ind_gps$outsidebuffer == "TRUE") &
        (ind_gps$count >= 12)
    ),
    "Excursion",
    "Not Excursion"
  )
  
}

##################################
#Step-selection Function Models

library(survival)

# create a blank list for model selection
mod.sel <- list()

# Global 3-way
mod.sel[[1]] <- clogit(#response variable
                       case ~ 
                       #movement covariates
                       steplength + log_sl + cos_ta +
                       #three-way interaction
                       Status*Sex*(Elevation + Distance.to.Herbaceous + Distance.to.Shrub + Distance.to.Water + Slope) +
                       #stratum
                       strata(step_id) + 
                       #cluster around individual
                       cluster(animal_id),
                       method = "approximate", data=all.data, model=TRUE)

# 2 way status only
mod.sel[[2]] <- clogit(#response variable
                       case ~ 
                       #movement covariates
                       steplength + log_sl + cos_ta +   
                       #two-way interaction with landscape familiarity status
                       Status*(Elevation + Distance.to.Herbaceous + Distance.to.Shrub + Distance.to.Water + Slope) +
                       #stratum
                       strata(step_id) + 
                       #cluster around individual
                       cluster(animal_id),
                       method = "approximate", data=all.data, model=TRUE)

# 2 way sex 
mod.sel[[3]] <- clogit(#response variable
                       case ~ 
                       #movement covariates
                       steplength + log_sl + cos_ta +
                       #two-way interaction with sex
                       Sex*(Elevation + Distance.to.Herbaceous + Distance.to.Shrub + Distance.to.Water + Slope) +
                       #stratum
                       strata(step_id) + 
                       #cluster around individual
                       cluster(animal_id),
                       method = "approximate", data=all.data, model=TRUE)


# null 
mod.sel[[4]] <- clogit(#response variable
                       case ~ 
                       #null intercept
                       1 + 
                       #stratum
                       strata(step_id),
                       #cluster around individual
                       cluster(animal_id),
                       method = "approximate", data = all.data, model = TRUE) 


# Name models
Model.names <- c("global", "global_subset_1", "global_subset_2", "null")

##################################

##################################
# Model selection based on AIC
##################################
library(AICcmodavg)
((out <- aictab(cand.set = mod.sel, modnames = Model.names, second.ord=F)))

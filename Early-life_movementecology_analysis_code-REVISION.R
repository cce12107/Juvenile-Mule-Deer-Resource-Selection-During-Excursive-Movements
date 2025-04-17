# Early-life movement ecology of a cervid: Implications for Chronic Wasting Disease Management Analysis Code
# Ellis et al. 2025 in Ecological Solutions & Evidence
# Calvin C. Ellis
# 02/18/2025

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
library(MuMIn)

### Males During Excursions ###
male.exc.global = clogit(#response variable
  case ~ 
    #movement covariates
    sl_ + log_sl_ + cos_ta_ +
    #landscape covariates
    Elevation + Distance.to.Herbaceous + Distance.to.Shrub + Distance.to.Water + Slope +
    #stratum
    strata(step_id_) + 
    #cluster around individual
    cluster(id),
  method = "approximate", data=male.excursion.df, model=TRUE)

#testing all combinations of models while including movement covariates in every iteration

options(na.action = 'na.fail')
male.exc.dredge <- dredge(male.exc.global,
                      subset = sl_ & log_sl_ & cos_ta_ & 
                        (Elevation | Distance.to.Herbaceous | Distance.to.Shrub | 
                           Distance.to.Water | Slope))

male.exc.model.results <- as.data.frame(male.exc.dredge)
male.exc.top = subset(male.exc.dredge, delta < 2) #extract competing models

### Males Within Range ###
male.natal.global = clogit(#response variable
  case ~ 
    #movement covariates
    sl_ + log_sl_ + cos_ta_ +
    #landscape covariates
    Elevation + Distance.to.Herbaceous + Distance.to.Shrub + Distance.to.Water + Slope +
    #stratum
    strata(step_id_) + 
    #cluster around individual
    cluster(id),
  method = "approximate", data=male.natal.df, model=TRUE)

#testing all combinations of models while including movement covariates in every iteration

options(na.action = 'na.fail')
male.natal.dredge <- dredge(male.natal.global,
                      subset = sl_ & log_sl_ & cos_ta_ & 
                        (Elevation | Distance.to.Herbaceous | Distance.to.Shrub | 
                           Distance.to.Water | Slope))

male.natal.model.results <- as.data.frame(male.natal.dredge)
male.natal.top = subset(male.natal.dredge, delta < 2) #extract competing models

### Females During Excursions ###
female.exc.global = clogit(#response variable
  case ~ 
    #movement covariates
    sl_ + log_sl_ + cos_ta_ +
    #landscape covariates
    Elevation + Distance.to.Herbaceous + Distance.to.Shrub + Distance.to.Water + Slope +
    #stratum
    strata(step_id_) + 
    #cluster around individual
    cluster(id),
  method = "approximate", data=female.excursion.df, model=TRUE)

#testing all combinations of models while including movement covariates in every iteration

options(na.action = 'na.fail')
female.exc.dredge <- dredge(female.exc.global,
                      subset = sl_ & log_sl_ & cos_ta_ & 
                        (Elevation | Distance.to.Herbaceous | Distance.to.Shrub | 
                           Distance.to.Water | Slope))

female.exc.model.results <- as.data.frame(female.exc.dredge)
female.exc.top = subset(female.exc.dredge, delta < 2) #extract competing models

### Females Within Range ###
female.natal.global = clogit(#response variable
  case ~ 
    #movement covariates
    sl_ + log_sl_ + cos_ta_ +
    #landscape covariates
    Elevation + Distance.to.Herbaceous + Distance.to.Shrub + Distance.to.Water + Slope +
    #stratum
    strata(step_id_) + 
    #cluster around individual
    cluster(id),
  method = "approximate", data=female.natal.df, model=TRUE)

#testing all combinations of models while including movement covariates in every iteration

options(na.action = 'na.fail')
female.natal.dredge <- dredge(female.natal.global,
                      subset = sl_ & log_sl_ & cos_ta_ & 
                        (Elevation | Distance.to.Herbaceous | Distance.to.Shrub | 
                           Distance.to.Water | Slope))

female.natal.model.results <- as.data.frame(female.natal.dredge)
female.natal.top = subset(female.natal.dredge, delta < 2) #extract competing models
##################################

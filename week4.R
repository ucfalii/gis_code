setwd("~/Desktop/GEOG0114/camden_theft")
library(sf)
library(nngeo)
install.packages(nngeo)
library(data.table)
library(tmap)
library(tidyverse)
install.packages(nngeo)
installed.packages(nngeo)
install.packages("nngeo")
library(nngeo)
camden_oas <- st_read('OAs_camden_2011.shp',crs=27700)
tm_shape(camden_oas)+tm_polygons()
tm_shape(camden_oas) +
  tm_borders(col='black') +
  tm_shape(camden_oas[camden_oas$OA11CD=='E00004174',]) +
  tm_fill(col='red')
chosen_oa <- 'E00004174'
chosen_oa_neighbours <- st_nn(st_geometry(st_centroid(camden_oas[camden_oas$OA11CD==chosen_oa,])), 
                              st_geometry(st_centroid(camden_oas)),
                              sparse = TRUE,
                              k = 50,
                              maxdist = 500) 
class(chosen_oa_neighbours)
neighbour_names <- camden_oas[chosen_oa_neighbours[[1]],]
neighbour_names <- neighbour_names$OA11CD
tm_shape(camden_oas) + 
  tm_borders() +
  # highlight only the neighbours
  tm_shape(camden_oas[camden_oas$OA11CD %in% neighbour_names,]) + 
  tm_fill(col = 'green') +
  # highlight only the chosen OA
  tm_shape(camden_oas[camden_oas$OA11CD==chosen_oa,]) + 
  tm_fill(col = 'red') +
  tm_shape(camden_oas) + 
  # overlay the borders
  tm_borders(col='black')
st_rook = function(a, b = a) st_relate(a, b, pattern = 'F***1****')
st_queen <- function(a, b = a) st_relate(a, b, pattern = 'F***T****')
chosen_oa_neighbours <- st_rook(st_geometry(camden_oas[camden_oas$OA11CD==chosen_oa,]), 
                                st_geometry(camden_oas))
neighbour_names <- camden_oas[chosen_oa_neighbours[[1]],]
neighbour_names <- neighbour_names$OA11CD
tm_shape(camden_oas) + 
  tm_borders() +
  # highlight only the neighbours
  tm_shape(camden_oas[camden_oas$OA11CD %in% neighbour_names,]) + 
  tm_fill(col = 'green') +
  tm_shape(camden_oas[camden_oas$OA11CD==chosen_oa,]) + 
  # highlight only the chosen OA
  tm_fill(col = 'red') +
  tm_shape(camden_oas) + 
  # overlay the borders
  tm_borders(col='black')

##Theft in Camden
camden_theft <- fread('crime/2019_camden_theft_from_person.csv')
camden_theft <- st_as_sf(camden_theft, coords = c('X','Y'), crs = 27700)
tm_shape(camden_oas) +
  tm_polygons() +
  tm_shape(camden_theft) +
  tm_dots()
# thefts in camden
camden_oas$n_thefts <- lengths(st_intersects(camden_oas, camden_theft))

# inspect
tm_shape(camden_oas) +
  tm_fill(col='n_thefts')

# square root of thefts
camden_oas$sqrt_n_thefts <- sqrt(camden_oas$n_thefts)

# inspect
tm_shape(camden_oas) +
  tm_fill(col='sqrt_n_thefts')

### Global Moran's 
library(sp)
library(spdep)
install.packages('spdep')
library(spdep)
class(camden_oas)
camden_oas_sp <- as_Spatial(camden_oas, IDs=camden_oas$OA11CD)
class(camden_oas_sp)
camden_oas_nb <- poly2nb(camden_oas_sp, row.names=camden_oas_sp$OA11CD)
class(camden_oas_nb)
str(camden_oas_nb,list.len = 10)
### create the list weights objected
nb_weights_list <- nb2listw(camden_oas_nb, style='W')
class(nb_weights_list)
### Moran's i
mi_value <- moran(camden_oas_sp$n_thefts,nb_weights_list,n=length(nb_weights_list$neighbours),S0=Szero(nb_weights_list))
mi_value
### run a Monte Carlo simulation 599 times
mc_model <- moran.mc(camden_oas_sp$n_thefts, nb_weights_list, nsim=599)
mc_model

### Local Moran's
# create an nb object
camden_oas_nb <- poly2nb(camden_oas_sp, row.names=camden_oas_sp$OA11CD)

# create the list weights object
nb_weights_list <- nb2listw(camden_oas_nb, style='W')

# Local Moran's I
local_moran_camden_oa_theft <- localmoran(camden_oas_sp$n_thefts, nb_weights_list)

# rescale
camden_oas_sp$scale_n_thefts <- scale(camden_oas_sp$n_thefts)

# create a spatial lag variable 
camden_oas_sp$lag_scale_n_thefts <- lag.listw(nb_weights_list, camden_oas_sp$scale_n_thefts)

# convert to sf
camden_oas_moran_stats <- st_as_sf(camden_oas_sp)

# set a significance value
sig_level <- 0.1

# classification with significance value
camden_oas_moran_stats$quad_sig <- ifelse(camden_oas_moran_stats$scale_n_thefts > 0 & 
                                            camden_oas_moran_stats$lag_scale_n_thefts > 0 & 
                                            local_moran_camden_oa_theft[,5] <= sig_level, 
                                          'high-high', 
                                          ifelse(camden_oas_moran_stats$scale_n_thefts <= 0 & 
                                                   camden_oas_moran_stats$lag_scale_n_thefts <= 0 & 
                                                   local_moran_camden_oa_theft[,5] <= sig_level, 
                                                 'low-low', 
                                                 ifelse(camden_oas_moran_stats$scale_n_thefts > 0 & 
                                                          camden_oas_moran_stats$lag_scale_n_thefts <= 0 & 
                                                          local_moran_camden_oa_theft[,5] <= sig_level, 
                                                        'high-low', 
                                                        ifelse(camden_oas_moran_stats$scale_n_thefts <= 0 & 
                                                                 camden_oas_moran_stats$lag_scale_n_thefts > 0 & 
                                                                 local_moran_camden_oa_theft[,5] <= sig_level, 
                                                               'low-high',
                                                               ifelse(local_moran_camden_oa_theft[,5] > sig_level, 
                                                                      'not-significant', 
                                                                      'not-significant')))))

# classification without significance value
camden_oas_moran_stats$quad_non_sig <- ifelse(camden_oas_moran_stats$scale_n_thefts > 0 & 
                                                camden_oas_moran_stats$lag_scale_n_thefts > 0, 
                                              'high-high', 
                                              ifelse(camden_oas_moran_stats$scale_n_thefts <= 0 & 
                                                       camden_oas_moran_stats$lag_scale_n_thefts <= 0, 
                                                     'low-low', 
                                                     ifelse(camden_oas_moran_stats$scale_n_thefts > 0 & 
                                                              camden_oas_moran_stats$lag_scale_n_thefts <= 0, 
                                                            'high-low', 
                                                            ifelse(camden_oas_moran_stats$scale_n_thefts <= 0 & 
                                                                     camden_oas_moran_stats$lag_scale_n_thefts > 0,
                                                                   'low-high',NA))))


# plot the results without the statistical significance
ggplot(camden_oas_moran_stats, aes(x = scale_n_thefts, 
                                   y = lag_scale_n_thefts, 
                                   color = quad_non_sig)) +
  geom_vline(xintercept = 0) + # plot vertical line
  geom_hline(yintercept = 0) + # plot horizontal line
  xlab('Scaled Thefts (n)') +
  ylab('Lagged Scaled Thefts (n)') +
  labs(colour='Relative to neighbours') +
  geom_point()



tm_shape(camden_oas)+tm_borders()+
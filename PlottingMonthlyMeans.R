
# Libraries ---------------------------------------------------------------
library(tidyverse)
library(raster)

# Load data ---------------------------------------------------------------
#Define file path
file_path <-'Temperature/Means/MonthlyMeans_2095-2100.nc'

#Load as raster brick
sst <- brick(file_path, stopIfNotEqualSpaced = F)


# Reproject data (Optional) -----------------------------------------------
#Reproject to South Polar Stereographic (EPSG:3976)
#Define projection
south_stereo <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

#Reproject raster brick
sst_reproj <- projectRaster(sst, crs =south_stereo)


# Plotting data -----------------------------------------------------------
#Base R option
plot(sst_reproj$X2100.01.16)
plot(sst$X2100.01.16)

#ggplot option - Requires data to be transformed into a data frame first
sst_df <- as.data.frame(sst$X2100.01.16, xy = T)
sst_df %>% 
  ggplot(aes(x, y, fill = X2100.01.16))+
  geom_raster()+
  scale_fill_gradient2(low = "#2c7bb6", mid = "#f7f7f7", 
                       high = "#ca0020", na.value = "#b6aca0")+
  theme_minimal()+
  labs(fill = "SST")

sst_df_reproj <- as.data.frame(sst_reproj$X2100.01.16, xy = T)
sst_df_reproj %>% 
  ggplot(aes(x, y, fill = X2100.01.16))+
  geom_raster()+
  scale_fill_gradient2(low = "#2c7bb6", mid = "#f7f7f7", 
                       high = "#ca0020", na.value = "#b6aca0")+
  theme_light()+
  labs(fill = "SST")

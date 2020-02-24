# Write Data from Stata-file to sf objects and export to shapefile readable with NetLogo GIS extension

library(sf)
library(tidyverse)
library(haven)
districts <- read_dta("TownData/ethnic_lsoa_town.2011_reduced.dta") %>% 
  mutate(town = as.character.factor(as_factor(town)), 
         lsoa11cd = as.character(lsoa11cd)) %>% 
  select(lsoa11cd, town, all, all1674valid, starts_with("whiteb_"), starts_with("asian_"), 
         starts_with("black_"), starts_with("othereth_"))
englandwales <- read_sf(dsn ="TownData/Lower_Layer_Super_Output_Areas_December_2011_Boundaries_EW_BSC/",
                        stringsAsFactors=FALSE) %>% 
  st_transform("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=m +no_defs") %>% 
  left_join(districts, by = c("LSOA11CD" = "lsoa11cd"))

for (TOWN in c("Bradford", "Leicester", "London", "Manchester", "Birmingham", "Blackburn", "Oldham","Wolverhampton","Slough","Walsall","Bristol", "Liverpool", "Leeds")) {
  englandwales %>% filter(town == TOWN) %>%
    st_write(dsn = paste0("shp_NetLogo/",TOWN), layer = TOWN, driver = "ESRI Shapefile", delete_dsn = TRUE)
}

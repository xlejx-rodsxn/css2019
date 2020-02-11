# Write Data from Stata-file to sf objects and export to shapefile readable with NetLogo GIS extension

library(sf)
library(haven)
districts <- read_dta("TownData/ethnic_lsoa_town.2011_reduced.dta") %>% 
  mutate(town = as.character.factor(as_factor(town)), 
         lsoa11cd = as.character(lsoa11cd)) %>% 
  select(lsoa11cd, town, all, all1674valid, starts_with("whiteb_"), starts_with("asian_"), 
         starts_with("black_"), starts_with("othereth_"))
for (town in c("Bradford", "Leicester", "London", "Manchester", "Birmingham")) {
  shp <- read_sf(dsn ="TownData/",layer = town, stringsAsFactors=FALSE) %>% 
    select(objectid, lsoa11cd, lsoa11nm, st_areasha, st_lengths, gor_name, town, geometry) %>% 
    left_join(districts, by = c("lsoa11cd", "town"))
  st_write(shp, dsn = paste0("shp_NetLogo/",town), layer = town, driver = "ESRI Shapefile", delete_dsn = TRUE)
}


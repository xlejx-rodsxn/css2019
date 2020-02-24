# Write Data from xlsx-File of statistik.bremen.de to sf objects and export to shapefile 
# readable with NetLogo GIS extension
library(sf)
library(tidyverse)
library(tidyxl)

bremen <- read_sf(dsn ="TownData/HB_pol_Grenzen_ETRS/ETRS/hb_ortsteile.shp",
                  stringsAsFactors=FALSE) %>% 
  st_transform("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=m +no_defs")

# Reading and wrangling xlsx Data 
cells <- xlsx_cells("TownData/Bremen_Datentabelle_Ortsteilatlas.xlsx") %>% 
  mutate(value = na_if(if_else(is_blank, "", if_else(data_type=="character",character,as.character(numeric))), "")) %>%
  select(sheet, row, col, value)
kennzahlen <- cells %>% filter(row==2, col>=3, !is.na(value)) %>% select(-row) %>% rename(Kennzahl=value)
jahre <- cells %>% filter(sheet == sheet, row == 3, col >= 3) %>% select(-row) %>% rename(Jahr = value) %>% 
  left_join(kennzahlen, by = c("sheet", "col")) %>% fill(Kennzahl)

# Compile full data in long format 
data_bremen <- cells %>% filter(row >= 4, col == 2) %>% select(-col) %>% 
  full_join(jahre, by = "sheet") %>% rename(Gebiet = value) %>% 
  left_join(cells, by = c("sheet", "row", "col")) %>% 
  select(Thema = sheet, Kennzahl, Jahr, Gebiet, value, -row, -col) %>% 
  mutate(value = as.numeric(value))

# Dataset on migration background
Migr <- data_bremen %>% select(-Thema) %>% 
  filter(Kennzahl == "Bevölkerungszahl insgesamt" |
           Kennzahl == "Anteil der Bevölkerung mit Migrationshintergrund an der Gesamtbevölkerung (%)" |
           Kennzahl == "Anteil der ausländischen Bevölkerung an der Gesamtbevölkerung (%)" |
           Kennzahl == "Anteil der Ausländer/-innen an der Bevölkerung mit Migationshintergrund (%)" |
           Kennzahl == "Anteil der (Spät-)Aussiedler/-innen an der Gesamtbevölkerung (%)" |
           Kennzahl == "Anteil der (Spät-)Aussiedler/-innen an der Bevölkerung mit Migationshintergrund (%)" |
           Kennzahl == "Anteil der Bevölkerung mit türkischem Migrationshintergrund an der Gesamtbevölkerung (%)" |
           Kennzahl == "Anteil der Bevölkerung mit türkischer Nationalität an der Gesamtbevölkerung (%)") %>% 
  pivot_wider(names_from = "Kennzahl", values_from = "value") %>% 
  select(Gebiet, Jahr, Pop = "Bevölkerungszahl insgesamt", 
         FracMigrBkgnd_Pop = "Anteil der Bevölkerung mit Migrationshintergrund an der Gesamtbevölkerung (%)", 
         FracForeign_Pop = "Anteil der ausländischen Bevölkerung an der Gesamtbevölkerung (%)", 
         FracForeign_MigrBkgnd = "Anteil der Ausländer/-innen an der Bevölkerung mit Migationshintergrund (%)",
         FracLateResettle_Pop = "Anteil der (Spät-)Aussiedler/-innen an der Gesamtbevölkerung (%)", 
         FracLateResettle_MigrBkgnd = "Anteil der (Spät-)Aussiedler/-innen an der Bevölkerung mit Migationshintergrund (%)", 
         FracTurkBkgnd_Pop = "Anteil der Bevölkerung mit türkischem Migrationshintergrund an der Gesamtbevölkerung (%)", 
         FracTurkNat_Pop = "Anteil der Bevölkerung mit türkischer Nationalität an der Gesamtbevölkerung (%)") %>% 
  arrange(Gebiet, Jahr) %>%  fill(Pop) %>% # fill gap years of biannual population statistics 
  mutate(MigrBkgnd = Pop * FracMigrBkgnd_Pop / 100, 
         Foreign = Pop * FracForeign_Pop / 100, 
         Foreign2 = MigrBkgnd * FracForeign_MigrBkgnd / 100, 
         LateResettle = Pop * FracLateResettle_Pop / 100, 
         LateResettle2 = MigrBkgnd * FracLateResettle_MigrBkgnd / 100, 
         TurkBkgnd = Pop * FracTurkBkgnd_Pop / 100, 
         TurkNat = Pop * FracTurkNat_Pop / 100,
         MigrBkgndOther = MigrBkgnd - TurkBkgnd - LateResettle) 

# Migrants Data for NetLogo
MigrNetlogo <- Migr %>% select(Gebiet, Jahr, Pop, LateResettle, TurkBkgnd, MigrBkgndOther) %>% 
  mutate(LateResettle = round(LateResettle), 
         TurkBkgnd = round(TurkBkgnd),
         MigrBkgndOther = round(MigrBkgndOther),
         Native = Pop - LateResettle - TurkBkgnd - MigrBkgndOther)
MigrNetlogo_shp <- bremen %>% left_join(filter(MigrNetlogo, Jahr == 2017), by = c("BEZ_GEM" = "Gebiet"))
MigrNetlogo_shp %>% filter(Pop > 1400 & !is.na(Pop) & 
  OBJID != "DENIAT010000jvyg" & OBJID != "DENIAT010000jvbn") %>% # Remove small shapes which are part of 2-shape region
  st_write(dsn = "shp_NetLogo/Bremen", layer = "Bremen", driver = "ESRI Shapefile", delete_dsn = TRUE)


MigrNetlogo <- Migr %>% select(Gebiet, Jahr, Pop, LateResettle, TurkBkgnd, MigrBkgndOther) %>% 
  mutate(LateResettle = round(LateResettle), 
         TurkBkgnd = round(TurkBkgnd),
         MigrBkgndOther = round(MigrBkgndOther),
         MigrBkgnd = MigrBkgndOther + TurkBkgnd + LateResettle, 
         Native = Pop - MigrBkgnd) %>% select(-LateResettle, -TurkBkgnd, -MigrBkgndOther)
MigrNetlogo_shp <- bremen %>% left_join(filter(MigrNetlogo, Jahr == 2017), by = c("BEZ_GEM" = "Gebiet"))
MigrNetlogo_shp %>% filter(Pop > 1400 & !is.na(Pop) & 
                             OBJID != "DENIAT010000jvyg" & OBJID != "DENIAT010000jvbn") %>% 
  st_write(dsn = "shp_NetLogo/Bremen", layer = "Bremen", driver = "ESRI Shapefile", delete_dsn = TRUE)

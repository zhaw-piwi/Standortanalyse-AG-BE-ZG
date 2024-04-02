

library(terra)
library(dplyr)
library(tidyr)
library(lubridate)
library(sf)
library(purrr)
library(glue)
# ```
# for years in 2012-2021:
#   for days in 01.04.year-30.09.year):
#   ((tas-10)+(tas_max-10))/2*1.045
# 
# ```

# TabsD: Daily mean temperature (1961 – present)
# TmaxD: Daily maximum temperature (1961 – present)

kantone <- read_sf("data/swissboundaries3d/swissBOUNDARIES3D_1_5_LV95_LN02.gpkg", "tlm_kantonsgebiet")[,"name"]

kantone_filter <- kantone |> 
  filter(name %in% c("Aargau", "Bern", "Zug"))

ncs <- list.files("data/Klimadaten_Feb24", "\\.nc", full.names = TRUE,recursive = TRUE)

files_df <- tibble(filename = ncs) |> 
  mutate(basename = basename(filename)) |> 
  tidyr::separate_wider_delim(cols = basename, "_", names = c("Variable",NA,"From","To")) |> 
  mutate(
    From = as.Date(From, format = "%Y%m%d"),
    To = as.Date(To, format = "%Y%m%d")
  )


files_recent <- files_df |> 
  mutate(year = year(From)) |> 
  filter(year >= 2012) |> 
  pivot_wider(names_from = Variable, values_from = filename)


pmap(files_recent, \(From, To, year, TabsD, TmaxD){
  # browser()
  tabsd <- rast(TabsD)
  tmaxd <- rast(TmaxD)
  
  tabsd_season <- tabsd[[month(time(tabsd)) >= 4 & month(time(tabsd)) <= 9]]
  tmaxd_season <- tmaxd[[month(time(tmaxd)) >= 4 & month(time(tmaxd)) <= 9]]
  
  K = 1.045
  HI_ch <- sum((tabsd_season+tmaxd_season)/2-10)*K
  
  HI <- mask(HI_ch, kantone_filter)
  
  writeRaster(HI, glue("data-out/HI-historic/{year}.tif"),overwrite = TRUE)
  
},.progress = TRUE)






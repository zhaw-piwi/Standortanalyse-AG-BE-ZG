

################################################################################
## Load libraries / Functions
################################################################################

library("terra")
library("dplyr")
library("tidyr")
library("lubridate")
library("sf")
library("purrr")
library("glue")
library("ggplot2")
library("readr")

terraOptions(progress=0) # default is 3
overwrite <- TRUE


filter_raster_year <- function(myrast, from, to){
  # both sides included. This is not R's default (see cut() and terra::classify()).
  # change at some point?
  library(lubridate)
  myrast[[(year(time(myrast)) >= from & year(time(myrast)) <= to)]]
}

filter_raster_month <- function(myrast, from, to){
  # both sides included. This is not R's default (see cut() and terra::classify()).
  # change at some point?
  library(lubridate)
  myrast[[month(time(myrast)) >= from & month(time(myrast)) <= to]]
}
################################################################################
## Load Data
################################################################################


kantone <- read_sf("data/swissboundaries3d/swissBOUNDARIES3D_1_5_LV95_LN02.gpkg", "tlm_kantonsgebiet")[,"name"]

kantone_filter <- kantone |> 
  filter(name %in% c("Aargau", "Bern", "Zug"))

ncs <- list.files("data/prognosedaten", "\\.nc", full.names = TRUE,recursive = TRUE)

files_df <- tibble(filename = ncs) |> 
  mutate(basename = basename(filename)) |> 
  separate_wider_delim(cols = basename, "_", names = c("Variable","GCM_RCP",NA,"Year",NA,NA,NA)) |> 
  separate_wider_delim(cols = Year, "-", names = c("From", "To")) |> 
  separate_wider_regex(cols = GCM_RCP, c(GCM = ".*", "-", RCP = ".*"))




################################################################################
## Calculate HI per year
################################################################################

# 2 Zeitpunkte: 2031-2040, 2051-2060 (muss nach dem import gefiltert werden)
# 3 Emissionsszenarien: RCP 2.6, 4.5 und 8.5 (alles szenarien)
# 3 Klimamodelle: DMI-HIRHAM-ECEARTH, MPICSC-REMO1-MPIESM, SMHI-RCA-HADGEM (alle modelle)

files_wide <- files_df |> 
  pivot_wider(names_from = Variable, values_from = filename)

files_wide$RCP |> unique()
files_wide$GCM |> unique()


# The full set of periods is 2014, 2031 and 2051. However this seems to crash R
# probably when trying to do all the raster calculations in momory. The next 
# time I need to do this, it would be better to write out intermediate files
From_sel <- c(2014, 2031,2051)
# From_sel <- c(2014) 

files_wide <- files_wide |> 
  expand_grid(From_sel) |> 
  mutate(To_sel = From_sel + 9) |> 
  select(-From, -To) |> 
  mutate(
    name = paste(GCM, RCP, From_sel, To_sel, sep = "_")
  )

kantone_filter_4326 <- st_transform(kantone_filter, 4326)

HI_lowres <- files_wide |> 
  pmap(\(GCM, RCP, tas, tasmax, From_sel, To_sel, name){
  # browser()
  
  tabsd <- rast(tas)
  tmaxd <- rast(tasmax)
  
  tabsd_years <- filter_raster_year(tabsd, From_sel, To_sel)
  tmaxd_years <- filter_raster_year(tmaxd, From_sel, To_sel)
    
  tabsd_season <- filter_raster_month(tmaxd_years, 4, 9)
  tmaxd_season <- filter_raster_month(tmaxd_years, 4, 9)
    
  K = 1.045
  years <- unique(year(time(tabsd_season)))
    
  HI_ch <- lapply(years, \(x){
    # browser()
    tabsd_x <- filter_raster_year(tabsd_season, x, x)
    tmaxd_x <- filter_raster_year(tmaxd_season, x, x)
    sum((tabsd_x+tmaxd_x)/2-10,na.rm = TRUE)*K
    }) |> 
      rast() |> 
      mean()
    
  # add some buffer to kantone_filter_4326 ?
  HI <- crop(HI_ch, kantone_filter_4326)
  
  HI2 <- terra::project(HI,"epsg:2056")
  
  names(HI2) <- name
  HI2
},.progress = TRUE)


# I'm creating a copy of the df, so I dont run into issues when running the 
# previous loop in an interactive session.
files_wide2 <- files_wide
files_wide2$HI_lowres <- HI_lowres



################################################################################
## Get Mean over all GSMs
################################################################################

# Getting the mean now is probably faster, with the downside that I don't have 
# the raw values per GCM

# check if all rasters have a reasonable value range
# sapply(files_wide2$HI_lowres, minmax)
# 
# files_wide3 <- files_wide2 |>
#   group_by(RCP, From_sel, To_sel) |>
#   group_map(\(x,y){
#     # browser()
#     y <- mutate(y, name = paste(RCP, From_sel, To_sel, sep = "_"))
#     z <- rast(x$HI_lowres) |> mean()
#     names(z) <- y$name
#     y$HI_lowres <- list(z)
#     y
#   }) |>
#   bind_rows()



################################################################################
## Downscale to 10m resolution via Elevation
################################################################################


swissaltiregio <- rast("data/swissAltiRegio/swissaltiregio_2056_5728.tif")

swissaltiregio_crop <- crop(swissaltiregio, files_wide2$HI_lowres[[1]])

swissaltiregio_res <- resample(swissaltiregio_crop, HI_lowres[[1]])

HI_downscale <- map(files_wide2$HI_lowres, \(x){
  browser()
  HI_norm <- x + swissaltiregio_res * 1.1895
  HI_down <- resample(HI_norm, swissaltiregio_crop)

  HI_up <- HI_down - swissaltiregio_crop * 1.1895
  HI_up
}, .progress = TRUE)

files_wide2$HI_downscale <- HI_downscale

################################################################################
## Correct HI via Exposition and Slope
################################################################################

swissaltiregio_slope_deg <- terrain(swissaltiregio_crop, v = "slope", unit = "degrees") |> 
  round() |> # reclassify with "is → becomes"
  as.int()

swissaltiregio_aspect <- terrain(swissaltiregio_crop, v = "aspect")

ausrichtung <- seq(0, 360, 45)

himmelsrichtungen <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")
himmelsrichtungen_int <- c(1:8,1)


classes_df <- tibble(himmelsrichtungen, himmelsrichtungen_int, ausrichtung) |> 
  mutate(
    von = ausrichtung - 45/2,
    von = ifelse(von < 0, 0, von),
    bis = ausrichtung + 45/2,
    bis = ifelse(bis > 360, 360, bis)
    ) 

classes <- classes_df |> 
  select(von, bis, himmelsrichtungen_int) |> 
  as.matrix()


swissaltiregio_aspect2 <- classify(swissaltiregio_aspect,classes, include.lowest = TRUE)

cls_df <- classes_df |> 
  select(himmelsrichtungen_int, himmelsrichtungen) |> 
  head(-1)# since N is twice

levels(swissaltiregio_aspect2) <- cls_df

korrekturfaktoren <- read_csv("korrekturfaktoren.csv") |> 
  mutate(Hangneigung = as.integer(Hangneigung))


# Since we need to reclassify by 2 conditions (slope and aspect), we create a
# new classification matrix that takes both factors into account (apect * 100 + slope)

classes_df2 <- classes_df |> 
  select(starts_with("himmelsrichtungen")) |> 
  head(-1)

korr_classify <- korrekturfaktoren |> 
  pivot_longer(-Hangneigung, names_to = "himmelsrichtungen", values_to = "faktor") |> 
  left_join(classes_df2, by = "himmelsrichtungen") |> 
  transmute(
    from = himmelsrichtungen_int*100 + Hangneigung,
    to = faktor
  )


swissaltiregio_aspect_slope <- swissaltiregio_aspect2*100+swissaltiregio_slope_deg

# This creates a "correction" raster based on the actual slope / aspect 
slope_aspect_korrektur <- classify(swissaltiregio_aspect_slope, korr_classify, others = NA)

# This correction raster is applied to a normalized version of the HI values
slope_aspect_korrektur_normalized <- (2000 * (1 - slope_aspect_korrektur))
HI_corrected <- map(files_wide2$HI_downscale, \(x){
  
  y <- x - slope_aspect_korrektur_normalized
  y[y < 0] <- 0
  y
}, .progress = TRUE)


files_wide2$HI_corrected <- HI_corrected


################################################################################
## Export
################################################################################


files_wide2 |> 
  select(HI_downscale, HI_corrected) |> 
  pmap(\(HI_downscale, HI_corrected){
    # browser()
  filename <- names(HI_downscale)
  HI_corrected <- as.int(HI_corrected)
  writeRaster(HI_corrected, filename = glue("data-out/HI-future/{filename}.tif"), overwrite = overwrite)
  
  # HI_downscale <- as.int(HI_downscale)
  # writeRaster(x, filename = glue("data-out/HI-future/uncorrected/{HI_downscale}.tif"), overwrite = overwrite)
}, .progress = TRUE)







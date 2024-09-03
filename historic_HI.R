

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

overwrite <- FALSE

terraOptions(progress = 0)

################################################################################
## Load Data
################################################################################


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
  filter(year >= 1971) |>
  pivot_wider(names_from = Variable, values_from = filename)


################################################################################
## Calculate HI per year
################################################################################



ret <- pmap(files_recent, \(From, To, year, TabsD, TmaxD){
  # browser()
  tabsd <- rast(TabsD)
  tmaxd <- rast(TmaxD)
  
  tabsd_season <- tabsd[[month(time(tabsd)) >= 4 & month(time(tabsd)) <= 9]]
  tmaxd_season <- tmaxd[[month(time(tmaxd)) >= 4 & month(time(tmaxd)) <= 9]]
  
  K = 1.045
  HI_ch <- sum((tabsd_season+tmaxd_season)/2-10)*K
  
  HI <- crop(HI_ch, kantone_filter)
  
  # writeRaster(HI_mueller_thurgau, glue("data-out/HI-historic/{year}.tif"),overwrite = overwrite)
  names(HI) <- year
  HI
},.progress = TRUE)


HI_all <- rast(ret)


################################################################################
## Mean HI per decade
################################################################################


from <- c(1974, 1984, 1994, 2004, 2014)

HI_av <- map(from, \(from){
  
  to <- from + 10
  
  years <- HI_all |> 
    names() |> 
    as.integer() 
  
  HI_mean <- mean(HI_all[[years >= from & years < to]])
  names(HI_mean) <- paste(from, to, sep = "-")
  HI_mean
}, .progress = TRUE)


################################################################################
## Downscale to 10m resolution via Elevation
################################################################################

# Das Downscaling wurde in zwei Schritten durchgeführt. Als erstes wurde eine Höhenkorrektur des HI
# durchgeführt. Dafür wurde das Ausgangsraster auf Meereshöhe normalisiert. Dabei wurde von einer
# Abnahme der Lufttemperatur von 6.5°C pro 100 Höhenmeter ausgegangen 10,11. Für den Huglin-Index
# entspricht das einer Abnahme von 1.1895 Punkten pro Höhenmeter (183 Tage * 0.0065 °C/m) und
# wurde wie folgt berechnet:

# HI-Raster normalisiert auf Meereshöhe (1km) = HI-Ausgangsraster (1km) + (Höhenlage (1km) * 1.1895)

# Danach wurde die räumliche Auflösung des normalisierten HI-Rasters mittels resampling von 1km auf
# 2m runterskaliert. Dabei wurden die resultierenden HI-Werte bilinear interpoliert. Im Anschluss wurde
# der HI-Wert für die tatsächliche Höhenlage in der 2m-Auflösung wie folgt berechnet:

# HI-Raster höhenkorrigiert (2m) = HI-Raster normalisiert auf Meereshöhe (2m) – (Höhenlage (2m) * 1.1895)

swissaltiregio <- rast("data/swissAltiRegio/swissaltiregio_2056_5728.tif")

swissaltiregio_crop <- crop(swissaltiregio, HI_av[[1]])

swissaltiregio_res <- resample(swissaltiregio_crop, HI_av[[1]])

HI_downscale <- map(HI_av, \(x){
  # browser()
  HI_norm <- x + swissaltiregio_res * 1.1895
  HI_down <- disagg(HI_norm, 100, method = "bilinear")
  HI_up <- HI_down - swissaltiregio_crop * 1.1895
  HI_up
}, .progress = TRUE)




################################################################################
## Correct HI via Exposition and Slope
################################################################################

swissaltiregio_slope_deg <- terrain(swissaltiregio_crop, v = "slope", unit = "degrees") |> 
  round() |> # dont reclassify with "from-to → becomes" but "is → becomes"
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
HI_corrected <- map(HI_downscale, \(x){
  
  y <- x - slope_aspect_korrektur_normalized
  y[y < 0] <- 0
  y
}, .progress = TRUE)


################################################################################
## Export
################################################################################



# Export the historic, downscaled HI values
map(HI_corrected, \(x){
  years <- names(x)
  if(overwrite){
    writeRaster(x, glue("data-out/HI-historic/corrected/{years}.tif"),overwrite = overwrite)
  }
}, .progress = TRUE)




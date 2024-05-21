

library(terra)
library(dplyr)
library(tidyr)
library(lubridate)
library(sf)
library(purrr)
library(glue)
library(ggplot2)

overwrite <- FALSE


# given a vector of length n, select k elements in a rolling /
# moving window fashion. This will produce a list of length
# n - k + 1
roll_select <- function(n, k){
  
  stopifnot(n >= k)
  lapply(seq_len(n-k+1), \(x){
    vec <- rep(FALSE,n)
    vec[x:(x+k-1)] <- TRUE
    vec
  })
}

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

# rolling mean
# HI_rollmean <- lapply(roll_select(nlyr(HI_all),11),\(x){
#   HI_sel <- HI_all[[x]]
#   year <- paste(range(as.numeric(names(HI_sel))),collapse = "-")
#   HI_mean <- mean(HI_sel)
#   names(HI_mean) <- year
#   HI_mean
# })

# HI_all |> 
#   names() |> 
#   as.integer() |> 
#   (\(x) x %% 10 == 3)()

from <- c(1974, 1984, 1994, 2004, 2014)

HI_av <- map(from, \(from){
  
  to <- from + 10
  
  years <- HI_all |> 
    names() |> 
    as.integer() 
  
  HI_mean <- mean(HI_all[[years >= from & years < to]])
  names(HI_mean) <- paste(from, to, sep = "-")
  HI_mean
})


swissaltiregio <- rast("data/swissAltiRegio/swissaltiregio_2056_5728.tif")

swissalti3d <- rast("data/swissAlti3D_mosaic.tif")
# dhm25 <- rast("data/DHM25_MM_ASCII_GRID/ASCII_GRID_1part/dhm25_grid_raster.asc")
# crs(dhm25) <- "epsg:21781"

# dhm25_2056 <- project(dhm25, "epsg:2056")
plot(swissaltiregio)

swissaltiregio_crop <- crop(swissaltiregio, HI_av[[1]])

swissaltiregio_res <- resample(swissaltiregio_crop, HI_av[[1]],"bilinear")

# swissaltiregio_res2 <- resample(swissaltiregio_crop, disagg(HI_av[[1]], 100), "bilinear")


HI_downscale <- map(HI_av, \(x){
  HI_norm <- x + swissaltiregio_res * 1.1895
  HI_down <- disagg(HI_norm, 100)
  HI_up <- HI_down - swissaltiregio_crop * 1.1895
  HI_up
}, .progress = TRUE)


# Get the full range of both rasters
lims <- range(minmax(HI_av[[1]]), minmax(HI_downscale[[1]]))

# plot both rasters with the range of both rasters
# (makes it easier to compare them)
plot(HI_av[[1]], range = lims)
plot(HI_downscale[[1]], range = lims)


# Export the historicm, downscaled HI values
map(HI_downscale, \(x){
  years <- names(x)
  if(overwrite) writeRaster(x, glue("data-out/HI-historic/{years}.tif"),overwrite = overwrite)
}, .progress = TRUE)


swissaltiregio_slope <- terrain(swissaltiregio_crop, v = "slope", unit = "radians")
swissaltiregio_slope_perc <- tan(swissaltiregio_slope) * 100

writeRaster(swissaltiregio_crop, "data-out/swissaltiregio.tif", overwrite = overwrite)
writeRaster(swissaltiregio_slope_perc, "data-out/slope.tif", overwrite = overwrite)


swissaltiregio_aspect <- terrain(swissaltiregio_crop, v = "aspect")
ausrichtung <- seq(0, 360, 45)

himmelsrichtungen <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")
names(ausrichtung) <- himmelsrichtungen
# this seems to be the simplest way:
himmelsrichtungen_int <- c(1:8,1)

classes <- matrix(
  c(
    ausrichtung-45/2,
    ausrichtung+45/2,
    himmelsrichtungen_int
  ),ncol = 3
)

swissaltiregio_aspect2 <- classify(swissaltiregio_aspect,classes)

cls_df <- tibble(id = himmelsrichtungen_int, azimut = himmelsrichtungen) |> 
  head(-1) # since N is twice

levels(swissaltiregio_aspect2) <- cls_df
writeRaster(swissaltiregio_aspect2, "data-out/exposition_class.tif", overwrite = overwrite)
writeRaster(swissaltiregio_aspect, "data-out/exposition.tif", overwrite = overwrite)



# 
# 
# 
# HI_sf2 <- do.call(rbind, HI_sf)
# 
# 
# 
# HI_sf3 <- HI_sf2 |> 
#   separate_wider_delim(year, "-",names = c("from","to"),cols_remove = FALSE) |> 
#   mutate(across(c(from,to), as.integer)) |> 
#   st_as_sf() |> 
#   group_by(from, to, year) |> 
#   summarise() |> 
#   arrange(desc(from))
# 
# write_sf(HI_sf3, "data-out/HI_sf3.gpkg")
# 
# # library(smoothr)
# # 
# # HI_chaikin <- smoothr::smooth(HI_sf3,method = "chaikin")
# # HI_chaikin_20 <- smoothr::smooth(HI_sf3,method = "chaikin", refinements = 10)
# # 
# # 
# # HI_ksmooth <- smoothr::smooth(HI_sf3,method = "ksmooth")
# # HI_spline <- smoothr::smooth(HI_sf3,method = "spline")
# # HI_densify <- smoothr::smooth(HI_sf3,method = "densify")
# 
# 
# ggplot(HI_sf3) +
#   geom_sf(data = kantone_filter, inherit.aes = FALSE) +
#   geom_sf(fill = "red", alpha = 0.3) +
#   # scale_fill_gradient(low = "red", high = "blue") +
#   facet_wrap(~year)
# 
# 
# 
# ggplot(HI_densify) +
#   geom_sf(data = kantone_filter, inherit.aes = FALSE) +
#   geom_sf(aes(fill = from, alpha = 0.3)) +
#   scale_fill_gradientn(colours  = RColorBrewer::brewer.pal(4, "Spectral"))  +
#   facet_wrap(~year)
# 

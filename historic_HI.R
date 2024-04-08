

library(terra)
library(dplyr)
library(tidyr)
library(lubridate)
library(sf)
library(purrr)
library(glue)
library(ggplot2)


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
  
  # writeRaster(HI_mueller_thurgau, glue("data-out/HI-historic/{year}.tif"),overwrite = TRUE)
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

HI_all |> 
  names() |> 
  as.integer() |> 
  (\(x) x %% 10 == 3)()

from <- c(1974, 1984, 1994, 2004, 2014)
to <- from + 9

HI_av <- map2(from, to, \(from, to){
  
  years <- HI_all |> 
    names() |> 
    as.integer() 
  
  HI_mean <- mean(HI_all[[years >= from & years <= to]])
  names(HI_mean) <- paste(from, to, sep = "-")
  HI_mean
})


swissalti3d <- rast("data/swissAlti3D_mosaic.tif")
dhm25 <- rast("data/DHM25_MM_ASCII_GRID/ASCII_GRID_1part/dhm25_grid_raster.asc")
crs(dhm25) <- "epsg:21781"

dhm25_2056 <- project(dhm25, "epsg:2056")
plot(swissalti3d)


dhm25_res <- resample(dhm25_2056, HI_av[[1]],"bilinear")

dhm25_res2 <- resample(dhm25_2056, disagg(HI_av[[1]], 40))


HI_downscale <- sapply(HI_av, \(x){
  HI_norm <- x+dhm25_res * 1.1895
  HI_down <- disagg(HI_norm, 40)
  HI_up <- HI_down - dhm25_res2 * 1.1895
  HI_up
})


lapply(HI_downscale, \(x){
  years <- names(x)
  writeRaster(x, glue("data-out/HI-historic/{years}.tif"))
})

writeRaster(dhm25_res2, "data-out/dhm25.tif")

dhm25_slope <- terrain(dhm25_res2, v = "slope", unit = "radians")
dhm25_slope_perc <- tan(dhm25_slope) * 100

dhm25_aspect <- terrain(dhm25_res2, v = "aspect")

writeRaster(dhm25_res2, "data-out/dhm25.tif", overwrite = TRUE)
writeRaster(dhm25_slope_perc, "data-out/slope.tif", overwrite = TRUE)

ausrichtung <- seq(0, 360, 45)

himmelsrichtungen <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")
names(ausrichtung) <- himmelsrichtungen

himmelsrichtungen_int <- c(1:8,1)

classes <- matrix(
  c(
    ausrichtung-45/2,
    ausrichtung+45/2,
    himmelsrichtungen_int
  ),ncol = 3
)

dhm25_aspect2 <- classify(dhm25_aspect,classes)

cls_df <- tibble(id = himmelsrichtungen_int, azimut = himmelsrichtungen) |> 
  head(-1)

levels(dhm25_aspect2) <- cls_df
writeRaster(dhm25_aspect2, "data-out/exposition_class.tif", overwrite = TRUE)
writeRaster(dhm25_aspect, "data-out/exposition.tif", overwrite = TRUE)




HI_sf <- lapply(HI_downscale, \(x){
  x
  
  HI_classify <-  terra::classify(x,c(-Inf, 1400, 2100,Inf))
  
  levs_df <- tibble(value = 0:2, category = c("zu tief", "optimal", "zu hoch"))
  
  levels(HI_classify) <- levs_df
  
  HI_sf <- st_as_sf(as.polygons(HI_classify))
  
  HI_optimal <- HI_sf |>
    filter(category == "optimal") |>
    mutate(year = names(x)) |> 
    dplyr::select(-category)
  
  HI_optimal
})

HI_sf2 <- do.call(rbind, HI_sf)



HI_sf3 <- HI_sf2 |> 
  separate_wider_delim(year, "-",names = c("from","to"),cols_remove = FALSE) |> 
  mutate(across(c(from,to), as.integer)) |> 
  st_as_sf() |> 
  group_by(from, to, year) |> 
  summarise() |> 
  arrange(desc(from))

write_sf(HI_sf3, "data-out/HI_sf3.gpkg")

# library(smoothr)
# 
# HI_chaikin <- smoothr::smooth(HI_sf3,method = "chaikin")
# HI_chaikin_20 <- smoothr::smooth(HI_sf3,method = "chaikin", refinements = 10)
# 
# 
# HI_ksmooth <- smoothr::smooth(HI_sf3,method = "ksmooth")
# HI_spline <- smoothr::smooth(HI_sf3,method = "spline")
# HI_densify <- smoothr::smooth(HI_sf3,method = "densify")


ggplot(HI_sf3) +
  geom_sf(data = kantone_filter, inherit.aes = FALSE) +
  geom_sf(fill = "red", alpha = 0.3) +
  # scale_fill_gradient(low = "red", high = "blue") +
  facet_wrap(~year)



ggplot(HI_densify) +
  geom_sf(data = kantone_filter, inherit.aes = FALSE) +
  geom_sf(aes(fill = from, alpha = 0.3)) +
  scale_fill_gradientn(colours  = RColorBrewer::brewer.pal(4, "Spectral"))  +
  facet_wrap(~year)


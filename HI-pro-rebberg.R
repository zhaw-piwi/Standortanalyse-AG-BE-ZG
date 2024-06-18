

library(sf)
library(dplyr)
library(tmap)
library(terra)
library(tidyr)
library(stringr)
library(ggplot2)
library(glue)

# lwb_rebbaukataster_lv95 = Kanton Bern und Zug
be_zg <- read_sf("data/lwb_rebbaukataster_lv95/shapefiles/lwb_rebbaukataster_v2_0/rebbaukataster.shp")
ag <- read_sf("data/AGIS.al_rebkataster__kanton_aargau__Shapefile/AGIS.al_rebkataster/al_rebkataster_20230629.shp")

rebberge <- rbind(
  select(be_zg),
  select(ag)
)

tmap_mode("view")


hi_historic <- list.files("data-out/HI-historic/corrected", "\\.tif", full.names = TRUE) |> 
  lapply(rast) |> 
  rast()

hi_future <- list.files("data-out/HI-future/corrected/", "\\.tif", full.names = TRUE) |> 
  lapply(rast) |> 
  rast()


extracted_hi_historic <- extract(hi_historic, rebberge, mean, na.rm = TRUE)
extracted_hi_future <- extract(hi_future, rebberge, mean, na.rm = TRUE)

extracted_hi_historic

HI_values <- full_join(extracted_hi_future, extracted_hi_historic, by = "ID") |> 
  pivot_longer(-ID, values_to = "HI") |> 
  tidyr::extract(name, regex = "(RCP\\d{2}_)*?(\\d{4})[_-](\\d{4})",into = c("RCP", "from", "to")) |> 
  mutate(
    RCP = str_remove(RCP, "_"),
    RCP = ifelse(RCP == "", "n.a.",RCP)
  ) |> 
  mutate(across(c(from, to, HI), as.integer))



count(HI_values, from, to)

HI_values |>
  filter(RCP %in% c("n.a.", "RCP45")) |>
  ggplot(aes(from, HI, group = ID)) +
  geom_line(alpha = .2) +
  scale_x_continuous(breaks = unique(HI_values$from),labels = \(x)glue("{x}-{x+9}")) +
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank())


tm_shape(rebberge) + tm_polygons()

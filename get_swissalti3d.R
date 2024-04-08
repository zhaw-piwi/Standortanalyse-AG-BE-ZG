library(readr)
library(tidyr)
library(sf)
library(dplyr)
library(purrr)
library(glue)
ullr2poly <- function(xmin, ymin, xmax, ymax){
  c(1,2,1,4,3,4,3,2,1,2) |>
    sapply(\(x)c(xmin, ymin, xmax, ymax)[x]) |>
    matrix(ncol = 2, byrow = TRUE) |>
    list() |>
    st_polygon()
}

urls <- read_csv("data/swissalti3d_2m_all.csv",col_names = c("URL"))


urls <- urls |> 
  mutate(basename = basename(URL)) |> 
  separate_wider_regex(basename,
                       patterns = c(
                         "swissalti3d_\\d{4}_",
                         E = "\\d{4}",
                         "-",
                         N = "\\d{4}",
                         "_\\d_2056_\\d{4}.tif"
                       )
                       )

urls <- urls |> 
  # slice(1:2) |>
  mutate(
    E = as.integer(E)*1000,
    N = as.integer(N)*1000,
    E2 = E + 1000,
    N2 = N + 1000
  )

urls$geom <- pmap(urls, \(URL, E, N, E2, N2){ullr2poly(E,N, E2, N2)})

urls <- urls |> 
  st_as_sf(crs = 2056)

kantone <- read_sf("data/swissboundaries3d/swissBOUNDARIES3D_1_5_LV95_LN02.gpkg", "tlm_kantonsgebiet")[,"name"]

kantone_filter <- kantone |> 
  filter(name %in% c("Aargau", "Bern", "Zug"))

bb <- st_as_sfc(st_bbox(kantone_filter))
centr = st_centroid(bb)
bb2 <- (bb - centr) * 1.1 + centr
st_crs(bb2) <- 2056

urls_sel <- urls[bb2,,]

map(urls_sel$URL, \(x){
  bn <- basename(x)
  destfile <- glue("data/swissAlti3D/{bn}")
  if(!file.exists(destfile)){
    tryCatch(
      download.file(x, destfile, mode = "wb",quiet = TRUE), 
      error = \(x){x}
    )  
  }
  
},.progress = TRUE)





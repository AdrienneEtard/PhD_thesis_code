## read shapefiles
X <- c("ggmap", "rgdal", "maptools", "rgeos", "dplyr", "raster", "maps", "geosphere", "sp", "sf", "lwgeom", "mapproj", "ggplot2")
lapply(X, library, character.only = TRUE)
rm(X)

Shpf <-  st_as_sf(readOGR(dsn="../../Data/Ecoregions/tnc_terr_ecoregions.shp"))
class(Shpf)
Realms <- Shpf["WWF_REALM2"]
plot(Realms)

Try <- Realms %>% 
  dplyr::group_by(WWF_REALM2) %>%
  summarise(Polygon=st_union(geometry)) %>%
  ungroup()

class(Try)
plot(Try[,])

## 
try <- st_read("../../Data/Ecoregions/tnc_terr_ecoregions.shp")


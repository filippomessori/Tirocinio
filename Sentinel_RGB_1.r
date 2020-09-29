##### Visuallizzare immagine Sentinel in RGB


setwd("/Users/fillo/Desktop/TIROCINIO")
library(raster)
library(rgdal)
library(RStoolbox)

# Seleziono le immagini per le bande a 10m, 20m, 60m : purtroppo non son riuscito a caricare le immagini su R nel formato ".jp2" quindi le ho convertite in formato ".jpg"

setwd("/Users/fillo/Desktop/TIROCINIO")


### Banda 10m:

sentinel_10m <- brick("T32TPS_20200905T101031_TCI_10m.jpg")

par(mfrow=c(1,2))
plotRGB(sentinel_10m, r=3, g=2, b=1, stretch="Lin")
plotRGB(sentinel_10m, r=4, g=3, b=2, stretch="Lin")


# Banda 20m:

sentinel_20m <- brick("T32TPS_20200905T101031_TCI_20m.jpg")

par(mfrow=c(1,2))
plotRGB(sentinel_20m, r=3, g=2, b=1, stretch="Lin")
plotRGB(sentinel_20m, r=4, g=3, b=2, stretch="Lin")


# Banda 60m:

sentinel_60m <- brick("T32TPS_20200905T101031_TCI_60m.jpg")

par(mfrow=c(1,2))
plotRGB(sentinel_60m, r=3, g=2, b=1, stretch="Lin")
plotRGB(sentinel_60m, r=4, g=3, b=2, stretch="Lin")



# Confronto le bande tra loro:
 
par(mfrow=c(3,2))
plotRGB(sentinel_10m, r=3, g=2, b=1, stretch="Lin")
plotRGB(sentinel_10m, r=4, g=3, b=2, stretch="Lin")
plotRGB(sentinel_20m, r=3, g=2, b=1, stretch="Lin")
plotRGB(sentinel_20m, r=4, g=3, b=2, stretch="Lin")
plotRGB(sentinel_60m, r=3, g=2, b=1, stretch="Lin")
plotRGB(sentinel_60m, r=4, g=3, b=2, stretch="Lin")



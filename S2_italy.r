## Applying Spectral Diversity Indexes on North Eastern Alps (IT) ##
## Daniele Da Re, Enrico Tordoni, Duccio Rocchini ##
## October 2020

# rm(list = ls())
Sys.setlocale("LC_ALL", "English")

#---- 1. Set wd, load packages and custom functions ----
setwd("F:/UCLouvain/Backup20191112/working_files/rasterdiv/MEE_tfe/S2_italy/")
# save.image("workspace_mee_tfe.RData")
# load("workspace_mee_tfe_DR.RData")

library(tidyverse)
library(rasterdiv)
library(rgeos)
library(RStoolbox)
library(RTools)
library(unix)
library(rgdal)
library(gdalUtils)
library(sf)
library(RColorBrewer)
library(rasterVis)
require("proxy")
require("raster")
require("svMisc")
require(pbapply)
require(pbmcapply)
library(sen2r)
library(rasterVis)
library(ggplot2)
library(gridExtra)

rescale <- function(x, x.min = NULL, x.max = NULL, new.min = NULL, new.max = NULL) {
  if(is.null(x.min)) x.min = min(x)
  if(is.null(x.max)) x.max = max(x)
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

#---- 2. Load Sentinel 2 bands and preprocess them ----
bbox=shapefile("alps_test_area.shp")
#we first need to convert the sentinel data from jp2 to tif
setwd("F:/UCLouvain/Backup20191112/working_files/rasterdiv/MEE_tfe/S2_italy/S2A_MSIL2A_20200905T101031_N0214_R022_T32TPS_20200905T130252/S2A_MSIL2A_20200905T101031_N0214_R022_T32TPS_20200905T130252.SAFE/GRANULE/L2A_T32TPS_A027189_20200905T101637/IMG_DATA/R10m/ ")

inpath="F:/UCLouvain/Backup20191112/working_files/rasterdiv/MEE_tfe/S2_italy/S2A_MSIL2A_20200905T101031_N0214_R022_T32TPS_20200905T130252/S2A_MSIL2A_20200905T101031_N0214_R022_T32TPS_20200905T130252.SAFE/GRANULE/L2A_T32TPS_A027189_20200905T101637/IMG_DATA/R10m/"
ls=list.files()

gdalcall = paste("gdalinfo",inpath, ls[3])
system2("C:/Program Files/QGIS 3.10/OSGeo4W.bat", args=gdalcall)

for(i in 1:length(ls)){
  # i=3
  INPUT_file=paste0(inpath, ls[i])
  outname=paste0(inpath, tools::file_path_sans_ext(ls[i]), ".tif")
  gdalcall=paste("gdal_translate", INPUT_file, outname,  "-co TILED=YES", "--config GDAL_CACHEMAX 1000")
  system2("C:/Program Files/QGIS 3.10/OSGeo4W.bat", args=gdalcall)
}

#now we load the .tif files and we crop them 
rlist=list.files(pattern=".tif", full.names = T)

s2_stack=stack(rlist[2:5])
s2_stack=crop(s2_stack, extent(bbox) )

names(s2_stack)<-c("blue", "green", "red", "nir")

plotRGB(s2_stack, r=3, g=2, b=1, stretch="Lin") #true colors
plotRGB(s2_stack, r=4, g=3, b=2, stretch="Lin") #false colors

setwd("F:/UCLouvain/Backup20191112/working_files/rasterdiv/MEE_tfe/S2_italy/")

#---- 3. Prepare NDVI ----
alps_ndvi= (s2_stack$nir-s2_stack$red)/(s2_stack$nir+s2_stack$red)
# alps_ndvi=raster("ndvi_alps_May9th2020.tif")
alps_ndvi=crop(alps_ndvi, bbox)
plot(alps_ndvi)

#rescale continuous ndvi to 8bit
alps_ndvi_recl= round(rescale(alps_ndvi, x.min = -1, x.max = 1, new.min = 0, new.max = 255),0)

alps_ndvi_recl= rescale(alps_ndvi, x.min = -1, x.max = 1, new.min = 0, new.max = 255)
alps_ndvi_recl
plot(alps_ndvi_recl)
storage.mode(alps_ndvi_recl[]) = "integer"

#---- 4. compute Spectral indexes on the whole image ----
r= alps_ndvi_recl

#Shannon's Diversity
sha <- Shannon(r,window=9,na.tolerance=0.9,np=3,cluster.type="SOCK")

#Rao's quadratic Entropy
rao <- Rao(r,window=9,na.tolerance=0.9,dist_m="euclidean",shannon=FALSE,np=3,cluster.type="SOCK")

#Renyi's Index
ren <- Renyi(r,window=9,alpha=seq(0,2),na.tolerance=0.9,np=3,cluster.type="SOCK")
ren_stack=stack(ren)
names_renyi=paste0("renyi_alpha_", seq(0,2))

#stack and save results
names_rstak=c(names(r), "Shannon", "Rao", names_renyi )

rstack_alps=stack(r, sha, rao, ren_stack)
names(rstack_alps)=names_rstak

plot(rstack_alps)

#---- 5. Plot outputs  ----
p0=ggRGB(s2_stack, r=3, g=2, b=1, stretch="Lin", 1, geom_raster = TRUE) +
  ggtitle("RGB") #true colors

p0.2=ggRGB(s2_stack, r=4, g=3, b=2, stretch="Lin")+
  ggtitle("RGB False colours") #true colors #false colors

p1 <- ggR(rstack_alps$ndvi_alps_May9th2020, 1, geom_raster = TRUE) +
  ggtitle("NDVI") +
  scale_fill_gradient(low='light green', high='dark green', na.value=NA)

#if you have a landcover/habitat shapefile 
# p2= ggplot(nc) + 
#   ggtitle("Land use") + geom_sf(aes(fill = newclass), lwd=0) + 
#   scale_fill_viridis_d(option='viridis')+
#   labs(fill="Classes")

p3 <- ggR(rstack_alps$Shannon, 1, geom_raster = TRUE, stretch="lin") +
  ggtitle("Shannon's H") + 
  scale_fill_gradient(low='yellow', high='blue', na.value=NA)

p4 <- ggR(rstack_alps$Rao, 1, geom_raster = TRUE, stretch="lin") + 
  ggtitle("Rao's Q")+
  scale_fill_gradient(low='yellow', high='blue', na.value=NA)

p5 <- ggR(rstack_alps$renyi_alpha_0, 1, geom_raster = TRUE, stretch="lin") + 
  ggtitle("Rényi (alpha=0)") + 
  scale_fill_gradient(low='yellow', high='blue', na.value=NA)

p9 <- ggR(rstack_alps$renyi_alpha_2, 1, geom_raster = TRUE, stretch="lin") + 
  ggtitle("Rényi (alpha=2)") + 
  scale_fill_gradient(low='yellow', high='blue', na.value=NA)

grid.arrange(p0, p0.2, p1, p3, p5, p9, p4, nrow = 2) # this needs griExtra

#plot rgb with ggplot 
# https://shekeine.github.io/visualization/2014/09/27/sfcc_rgb_in_R

#---- 6. Postprocessing indexes ----
ncogr_m=shapefile("YourLandCoverShapefile.shp")
df_extr=raster::extract(rstack_alps, ncogr_m)
names(df_extr)=ncogr_m$class

box_df=data.frame()
for(i in 1:length(df_extr)){
  # i=1
  df_tmp=df_extr[[i]]
  df_tmp=df_tmp%>%
    as_tibble() %>% 
    rename("Shannon's H"=layer.1, "Rao's Q"=layer.2,  "Rényi (alpha=0)"=layer.3,  "Rényi (alpha=2)"=layer.4) %>% 
    mutate(class= ncogr_m$class[i]) %>%
    pivot_longer(!class)
  box_df=rbind(box_df, df_tmp)
}

box_df %>% 
  filter(name!="Rao's Q") %>% 
  ggplot(aes(class, value, fill=name))+
  geom_boxplot() +
  # geom_violin()+
  xlab("Land Use classes") + ylab("Value")+
  labs(fill="Spectral diversity indexes")+
  theme_light()+
  theme(text = element_text(size=12),legend.text = element_text(face = "italic"),  axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), legend.position = "bottom")

#for plotting
# outname=paste0("imgs/", "Agri_prop_ratio", Sys.Date(), ".png")
# png(outname,width = 20, height = 15, units = "cm", res=300)




####################################################### FILIPPO MESSORI 

setwd("/Users/fillo/Desktop/TIROCINIO/S2_italy/S2A_MSIL2A_20200905T101031_N0214_R022_T32TPS_20200905T130252/S2A_MSIL2A_20200905T101031_N0214_R022_T32TPS_20200905T130252.SAFE/GRANULE/L2A_T32TPS_A027189_20200905T101637/IMG_DATA/R10m/ ")


# Faccio un crop per analizzare RGB, NDVI e indici di una zona vicino al ghiacciaio




# RGB:

extension_glacier <- c(630000, 640000, 5170000, 5180000)

s2_glacier_stack=crop(s2_stack, extension_glacier)

plotRGB(s2_glacier_stack, r=3, g=2, b=1, stretch="Lin") #true colors

plotRGB(s2_glacier_stack, r=4, g=3, b=2, stretch="Lin") # false colors


# NDVI:

alps_glaciers_ndvi=crop(alps_ndvi, extension_glacier)

plot(alps_glaciers_ndvi)

alps_glaciers_ndvi_recl= round(rescale(alps_glaciers_ndvi, x.min = -1, x.max = 1, new.min = 0, new.max = 255),0)

alps_glaciers_ndvi_recl= rescale(alps_glaciers_ndvi, x.min = -1, x.max = 1, new.min = 0, new.max = 255)

plot(alps_glaciers_ndvi_recl)

storage.mode(alps_glaciers_ndvi_recl[]) = "integer"

g <- alps_glaciers_ndvi_recl


#### INDEX:

# Shannon 
sha <- Shannon(g,window=9,na.tolerance=0.9,np=3,cluster.type="SOCK")

# RAO
rao_glaciers <- Rao(g,window=9,na.tolerance=0.9,dist_m="euclidean",shannon=FALSE,np=3,cluster.type="SOCK")

#Renyi
ren_glaciers <- Renyi(g,window=9,alpha=seq(0,2),na.tolerance=0.9,np=3,cluster.type="SOCK")
ren_stack=stack(ren)
names_renyi=paste0("renyi_alpha_", seq(0,2))


# Plot Index
names_rstak=c(names(r), "Shannon", "Rao", names_renyi )
rstack_alps=stack(r, sha, rao, ren_stack)
names(rstack_alps)=names_rstak
rstack_alps_glacier=crop(rstack_alps, extension_glacier) # Crop zona interesse
plot(rstack_alps_glacier)


# Analisi e Plot finale: RGB, NDVI e Indici


p0_glacier=ggRGB(s2_glacier_stack, r=3, g=2, b=1, stretch="Lin", 1, geom_raster = TRUE) + ggtitle("RGB glacier")

p02_glacier=ggRGB(s2_glacier_stack, r=4, g=3, b=2, stretch="Lin", 1, geom_raster = TRUE) + ggtitle("RGB glacier false colours")

p1_glacier <- ggR(rstack_alps_glacier$layer, 1, geom_raster = TRUE) + ggtitle("NDVI glacier") + scale_fill_gradient(low='light green', high='dark green', na.value=NA)
# NDVI

p3_glacier <- ggR(rstack_alps_glacier$Shannon, 1, geom_raster = TRUE, stretch="lin") + ggtitle("Shannon's H glaciers") + scale_fill_gradient(low='yellow', high='blue', na.value=NA)


p4_glacier <- ggR(rstack_alps_glacier$Rao, 1, geom_raster = TRUE, stretch="lin") +  ggtitle("Rao's Q glacier")+ scale_fill_gradient(low='yellow', high='blue', na.value=NA)

p5_glacier <- ggR(rstack_alps_glacier$renyi_alpha_0, 1, geom_raster = TRUE, stretch="lin") + ggtitle("Rényi (alpha=0) glacier") + scale_fill_gradient(low='yellow', high='blue', na.value=NA)

p9_glacier <- ggR(rstack_alps_glacier$renyi_alpha_2, 1, geom_raster = TRUE, stretch="lin") + ggtitle("Rényi (alpha=2) glacier") + scale_fill_gradient(low='yellow', high='blue', na.value=NA)

grid.arrange(p0_glacier, p02_glacier, p1_glacier, p3_glacier, p5_glacier, p9_glacier,p4_glacier, nrow= 2)


####

# Faccio un crop per analizzare RGB, NDVI e indici di una zona a valle ( osservare diversità vicino centro abitato e agricoltura )

extension_valle <- c(635000, 650000, 5155000, 5170000)
s2_valle_stack=crop(s2_stack, extension_valle)

# RGB valle:

plotRGB(s2_valle_stack, r=3, g=2, b=1, stretch="Lin") # true colors
plotRGB(s2_valle_stack, r=4, g=3, b=2, stretch="Lin") # false colors



# NDVI valle:

alps_valle_ndvi=crop(alps_ndvi, extension_valle)
plot(alps_valle_ndvi)
alps_valle_ndvi_recl= round(rescale(alps_valle_ndvi, x.min = -1, x.max = 1, new.min = 0, new.max = 255),0)
alps_valle_ndvi_recl= rescale(alps_valle_ndvi, x.min = -1, x.max = 1, new.min = 0, new.max = 255)
plot(alps_valle_ndvi_recl)
storage.mode(alps_valle_ndvi_recl[]) = "integer"
v <- alps_valle_ndvi_recl


### INDEX: 
# Recupero lo stack fatto in prescedenza "rstack_alps=stack(r, sha, rao, ren_stack)", e faccio un crop sulla nuova estensione

rstack_alps_valle=crop(rstack_alps, extension_valle)
plot(rstack_alps_valle)


# Analisi e Plot finale: RGB, NDVI e Indici

p0_valle=ggRGB(s2_valle_stack, r=3, g=2, b=1, stretch="Lin", 1, geom_raster = TRUE) + ggtitle("RGB valle")

p02_valle=ggRGB(s2_valle_stack, r=4, g=3, b=2, stretch="Lin", 1, geom_raster = TRUE) + ggtitle("RGB valle false colours")

p1_valle <- ggR(rstack_alps_valle$layer, 1, geom_raster = TRUE) + ggtitle("NDVI valle") + scale_fill_gradient(low='light green', high='dark green', na.value=NA)

p3_valle <- ggR(rstack_alps_valle$Shannon, 1, geom_raster = TRUE, stretch="lin") + ggtitle("Shannon's H valle") + scale_fill_gradient(low='yellow', high='blue', na.value=NA)

p4_valle <- ggR(rstack_alps_valle$Rao, 1, geom_raster = TRUE, stretch="lin") +  ggtitle("Rao's Q valle")+ scale_fill_gradient(low='yellow', high='blue', na.value=NA)

p5_valle <- ggR(rstack_alps_valle$renyi_alpha_0, 1, geom_raster = TRUE, stretch="lin") + ggtitle("Rényi (alpha=0) valle") + scale_fill_gradient(low='yellow', high='blue', na.value=NA)

p9_valle <- ggR(rstack_alps_valle$renyi_alpha_2, 1, geom_raster = TRUE, stretch="lin") + ggtitle("Rényi (alpha=2) valle") + scale_fill_gradient(low='yellow', high='blue', na.value=NA)

grid.arrange(p0_valle, p02_valle, p1_valle, p3_valle, p5_valle, p9_valle,p4_valle, nrow= 2)

#####


# Adesso confronto tra loro i singoli RGB, NDVI e indici:


grid.arrange(p0, p0.2, p0_glacier, p02_glacier, p0_valle, p02_valle, nrow= 3) # RGB

grid.arrange(p1, p1_glacier, p1_valle, nrow= 1)  # NDVI

grid.arrange(p3, p3_glacier, p3_valle, p4, p4_glacier, p4_valle, p5, p5_glacier, p5_valle, p9, p9_glacier, p9_valle, nrow= 4) # Index



# Aggiungo Renyi con alpha = 1

ren_stack=stack(ren)
names_renyi=paste0("renyi_alpha_", seq(0,1,2))
names_rstak=c(names(r), "Shannon", "Rao", names_renyi )
rstack_alps=stack(r, sha, rao, ren_stack)

p6 <- ggR(rstack_alps$renyi_alpha_1, 1, geom_raster = TRUE, stretch="lin") + ggtitle("Rényi (alpha=1)") + scale_fill_gradient(low='yellow', high='blue', na.value=NA)

p6_valle <- ggR(rstack_alps_valle$renyi_alpha_1, 1, geom_raster = TRUE, stretch="lin") + ggtitle("Rényi (alpha=1) valle") + scale_fill_gradient(low='yellow', high='blue', na.value=NA)

p6_glacier <- ggR(rstack_alps_glacier$renyi_alpha_1, 1, geom_raster = TRUE, stretch="lin") + ggtitle("Rényi (alpha=1) glacier") + scale_fill_gradient(low='yellow', high='blue', na.value=NA)


# Total area (update)

grid.arrange(p0, p0.2, p1, p3, p5, p6, p9, p4, nrow = 2)


# Glacier area (update)


grid.arrange(p0_glacier, p02_glacier, p1_glacier, p3_glacier, p5_glacier, p6_glacier, p9_glacier,p4_glacier, nrow= 2)


# Valle (update)

grid.arrange(p0_valle, p02_valle, p1_valle, p3_valle, p5_valle, p6_valle, p9_glacier, p4_valle, nrow= 2)


# INDEX (update)

grid.arrange(p3, p3_glacier, p3_valle, p4, p4_glacier, p4_valle, p5, p5_glacier, p5_valle, p6, p6_glacier, p6_valle, p9, p9_glacier, p9_valle, nrow= 5)



####### Renyi 1-10
setwd("~/Desktop/TIROCINIO/R/S2_italy/S2A_MSIL2A_20200905T101031_N0214_R022_T32TPS_20200905T130252.SAFE/GRANULE/L2A_T32TPS_A027189_20200905T101637/IMG_DATA/sent_renyi_duccio")


ren_tot <- lapply(list.files(pattern="tif"), function(x) raster(x)/1000)
pal <- colorRampPalette(c("purple","blue","cyan","green","yellow","red"), bias=3)
cuts <- seq(0,5.017,length=20)
renyi <- stack(ren_tot)
r.range <- c(minValue(renyi), maxValue(renyi))

png("~/renyi_indx_sentinel.png",width = 480*4, height = 480*2,pointsize=20)
plot(stack(ren_tot),col=pal(length(cuts)),breaks=cuts, legend=FALSE, axes=FALSE)
plot(renyi[[1]], legend.only=TRUE, horizontal=FALSE, col=pal(length(cuts)), legend.width=2, legend.shrink=0.5, axis.args=list(at=seq(0,5,1), labels=seq(0, 5, 1), cex.axis=2), legend.args=list(text="Renyi's Index", side=3, font=2, line=0.5, cex=2),smallplot=c(.82,.85, .03,.3)); par(mar = par("mar"))
dev.off()

# bias più alto di 1 per dare più risoluzione ai colori alla fine della paletta, minore di 1 per più risoluzione all'inizio. 
# cuts è il numero di classi della paletta

# /1000 dato che ho moltiplicato (*1000) e trasformato in interi i raster per diminuirne le dimensioni .


par(mfrow=c(3,4))
                  
plotRGB(s2_stack, r=3, g=2, b=1, stretch="Lin", axes=TRUE, main= "RGB", xlab= "x", ylab= "y") # RGB true colors
plotRGB(s2_stack, r=4, g=3, b=2, stretch="Lin", axes=TRUE, main= "RGB false colours", xlab= "x", ylab= "y") # RGB false colors
plot(alps_ndvi_recl, col= col_ndvi, xlab="x", ylab="y", main= "NDVI") # NDVI
plot(rstack_alps$Shannon, col=pal(length(cuts)),breaks=cuts, legend=T, axes=T, xlab="x", ylab= "y", main= "Shannon")
                  

                  
plot(renyi$Renyi_alpha_0, col=pal(length(cuts)),breaks=cuts, legend=T, axes=T, xlab="x", ylab= "y", main= "Renyi (alpha=0)") # alpha=0
plot(renyi$Renyi_alpha_1, col=pal(length(cuts)),breaks=cuts, legend=T, axes=T, xlab="x", ylab= "y", main= "Renyi (alpha=1)")# alpha=1
plot(renyi$Renyi_alpha_2, col=pal(length(cuts)),breaks=cuts, legend=T, axes=T, xlab="x", ylab= "y", main= "Renyi (alpha=2)")# alpha=2
plot(renyi$Renyi_alpha_3, col=pal(length(cuts)),breaks=cuts, legend=T, axes=T, xlab="x", ylab= "y", main= "Renyi (alpha=3)")# alpha=3
 
                  
plot(renyi$Renyi_alpha_4, col=pal(length(cuts)),breaks=cuts, legend=T, axes=T, xlab="x", ylab= "y", main= "Renyi (alpha=4)") # alpha=4
plot(renyi$Renyi_alpha_5, col=pal(length(cuts)),breaks=cuts, legend=T, axes=T, xlab="x", ylab= "y", main= "Renyi (alpha=5)") # alpha=5
plot(renyi$Renyi_alpha_10, col=pal(length(cuts)),breaks=cuts, legend=T, axes=T, xlab="x", ylab= "y", main= "Renyi (alpha=10)") # alpha=10
plot(rstack_alps$Rao, col=pal(length(cuts)),breaks=cuts, legend=T, axes=T, xlab="x", ylab= "y", main= "Rao") # Rao
                  
                  
                  

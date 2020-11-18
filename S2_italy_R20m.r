setwd("~/Desktop/TIROCINIO/R/S2_italy/S2A_MSIL2A_20200905T101031_N0214_R022_T32TPS_20200905T130252.SAFE/GRANULE/L2A_T32TPS_A027189_20200905T101637/IMG_DATA/R20m/R20m_update")

# RGB
rlist=list.files(pattern=".tif", full.names = T)
stack_20m <- stack(rlist[1:4])
names(stack_20m)<-c("blue", "green", "red", "nir")
extension_sgiustina <- c(645000, 665000, 5130000, 5145000)
stack_20m_sgiustina <- crop(stack_20m, extension_sgiustina)
plotRGB(stack_20m_sgiustina, r=4, g=3, b=2, stretch="Lin")
p0=ggRGB(stack_20m_sgiustina, r=3, g=2, b=1, stretch="Lin", 1, geom_raster = TRUE) + ggtitle("RGB Lago Santa Giustina")
p0.2=ggRGB(stack_20m_sgiustina, r=4, g=3, b=2, stretch="Lin", 1, geom_raster = TRUE) + ggtitle("RGB Lago false colours Santa Giustina")


# NDVI

ndvi_20m_sgiustina <- (stack_20m_sgiustina$nir - stack_20m_sgiustina$red)/(stack_20m_sgiustina$nir + stack_20m_sgiustina$red)
ndvi_20m_sgiustina_recl<- round(rescale(ndvi_20m_sgiustina, x.min = -1, x.max = 1, new.min = 0, new.max = 255),0)
ndvi_20m_sgiustina_recl= rescale(ndvi_20m_sgiustina, x.min = -1, x.max = 1, new.min = 0, new.max = 255)
plot(ndvi_20m_sgiustina_recl)

r = ndvi_20m_sgiustina_recl

# Index

sha <- Shannon(r,window=9,na.tolerance=0.9,np=3,cluster.type="SOCK")
rao <- Rao(r,window=9,na.tolerance=0.9,dist_m="euclidean",shannon=FALSE,np=3,cluster.type="SOCK")
ren <- Renyi(r,window=9,alpha=seq(0,2),na.tolerance=0.9,np=3,cluster.type="SOCK")
ren_stack=stack(ren)
names_renyi=paste0("renyi_alpha_", seq(0,1,2))




p1 <- ggR(ndvi_20m_sgiustina_recl, 1, geom_raster = TRUE) + ggtitle("NDVI") + scale_fill_gradient(low='light green', high='dark green', na.value=NA)
p3 <- ggR(sha, 1, geom_raster = TRUE, stretch="none") + ggtitle("Shannon's H") + scale_fill_gradient(low='yellow', high='blue', na.value=NA)
p4 <- ggR(rao, 1, geom_raster = TRUE, stretch="none") + ggtitle("Rao's Q") + scale_fill_gradient(low='yellow', high='blue', na.value=NA)
p5 <- ggR(ren_stack$Renyi_alpha_0, 1, geom_raster = TRUE, stretch="none") + ggtitle("Rényi (alpha=0)") + scale_fill_gradient(low='yellow', high='blue', na.value=NA)
p6 <- ggR(ren_stack$Shannon_Renyi_alpha_1, 1, geom_raster = TRUE, stretch="none") + ggtitle("Rényi (alpha=1)") + scale_fill_gradient(low='yellow', high='blue', na.value=NA)
p9 <- ggR(ren_stack$Renyi_alpha_2, 1, geom_raster = TRUE, stretch="none") + ggtitle("Rényi (alpha=2)") + scale_fill_gradient(low='yellow', high='blue', na.value=NA)



grid.arrange(p0, p0.2, p1, p3, p5, p6, p9, p4, nrow = 2) # Grafico

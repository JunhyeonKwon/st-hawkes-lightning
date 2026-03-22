########################################
##Packages
########################################
library(spacetime)
library(gstat)
#unloadNamespace(c('raster', 'terra')) 
require(xts)
require(forecast)
#for maps
require(maps)    # Provides functions that let us plot the maps
require(mapdata)    # Contains the hi-resolution points that mark out the countries.
require(sp)
require(rgdal)
require(maptools)
#require(rgeos)
#require(GeoXp)
require(ggplot2) #to use fortify function
library(shapefiles) #to simplyfy shape files
##For kriging
require(geoR) #for kriging
require(rrcov)
require(automap)

library(tidyverse)
library(animation)

library(lubridate)
library(sf)
library(spatstat)

library(PtProcess)
library(misc3d) #for kde3d fct
library(scatterplot3d)

library("chron")
library("fields")
#library("doMC")
library(ggplot2)
library(ggmap)
########################################
##Read data
########################################
library(xts)
library(lubridate)
#Lightning <- readRDS("~/Library/CloudStorage/Dropbox/Data/Lightning/Lightningevent.RDS")
Lightning <- readRDS("~/Library/CloudStorage/Dropbox/Data/Lightning/Lightningevent_2022.RDS")
#Lightning <- readRDS("C:\\Users\\HYU\\Dropbox\\Data\\Lightning\\Lightningevent.RDS")
unique_date <- unique(Lightning$data$일시)
#plot(Lightning$data$경도[which(Lightning$data$일시==unique_date[3])], Lightning$data$위도[which(Lightning$data$일시==unique_date[3])])

w0 <- Lightning$chull
result_mat3 <- Lightning$info

m <- length(unique(Lightning$data$event))
Lightning_data_list <- list()

for(j in 1:m){
  Lightning_data <- Lightning$data[which(Lightning$data$event==j),]
  time_spans <- unique(Lightning_data$일시)
  time_span <- c(time_spans[1], time_spans[length(time_spans)])
  t_span <- 6*(length(time_spans)-1)-2
  Lightning_data_df <- data.frame(long=double(), lat=double(), time=POSIXct(), tcount=double())
  for(k in 1:t_span){
    start_time <- time_span[1] + (k-1)*600
    end_time <- time_span[1] + (k-1)*600 + 1200
    which_index <- which(Lightning_data$시분초>=start_time & Lightning_data$시분초<=end_time)
    if(length(which_index)>0){
      Lightning_data_df <- rbind(Lightning_data_df,
                                 data.frame(long=mean(Lightning_data$경도[which_index]),
                                            lat=mean(Lightning_data$위도[which_index]),
                                            time=start_time+600,
                                            tcount=k))
    }
  }
  Lightning_data_list[[j]] <- Lightning_data_df
}

library(princurve)
par(mfrow=c(6,6))
par(mar=c(1.1,1.1,1.1,1.1))
for(i in 1:length(Lightning_data_list)){
  fit<-principal_curve(data.matrix(Lightning_data_list[[i]][,c("long", "lat", "tcount")]))
  #plot(fit,xlab="",ylab="")
  col_pal <- rev(heat.colors(length(Lightning_data_list[[i]]$long)))
  data_col <- col_pal[Lightning_data_list[[i]]$tcount]
  plot(Lightning_data_list[[i]]$long, Lightning_data_list[[i]]$lat, pch=16, col=data_col, xlab="", ylab="")
  lines(fit)
}


library(ggplot2)
library(ggmap)
library(patchwork) # install.packages("patchwork")
#NEED google keys
place <- "cheongju"
google <- get_googlemap(place, zoom = 10)
ggmap(google)

for(i in 1:34){
  if(i==1){
    Lightning_data_list2 <- Lightning_data_list[[i]]
  }else{
    Lightning_data_list2 <- rbind(Lightning_data_list2, Lightning_data_list[[i]])
  }
}

########################################
##Whole plot
########################################
bbox3 <- bb2bbox(attr(google, "bb"))
bbox3[1] <- min(Lightning_data_list2$long) - 1
bbox3[2] <- min(Lightning_data_list2$lat) - 1
bbox3[3] <- max(Lightning_data_list2$long) + 1
bbox3[4] <- max(Lightning_data_list2$lat) + 1

map_list <- lapply(1:34, function(i) {
  fit<-principal_curve(data.matrix(Lightning_data_list[[i]][,c("long", "lat", "tcount")]))
  curve_df <- as.data.frame(fit$s)
  colnames(curve_df) <- c("x", "y", "t")
  #curve_df <- curve_df[order(fit$lambda), ]
  
  map <-get_stadiamap(bbox3, maptype = "stamen_terrain_background", zoom=5) 
  ggmap(map)   + geom_point(data=Lightning_data_list[[i]], aes(x=long,y=lat, color=tcount), size=0.5) + geom_path(data = curve_df, aes(x = x, y = y), color="black", size=1) + scale_color_gradientn(colours = rainbow(20)) + labs(title=paste0(i), x = "", y = "")  + theme(legend.position = "none", axis.text = element_blank()) 
  #+ geom_sf(data=w0, inherit.aes = FALSE, fill=NA, color="black")
})


final_plot <- wrap_plots(map_list, ncol = 6)
#ggsave(paste0("PC-",0,".pdf"), plot=final_plot, width=10, height=10)

########################################
##Individual plot
########################################
bbox <- bb2bbox(attr(google, "bb"))
bbox[1] <- min(w0$geometry[[1]][[1]][,1]) - 0.25
bbox[2] <- min(w0$geometry[[1]][[1]][,2]) - 0.25
bbox[3] <- max(w0$geometry[[1]][[1]][,1]) + 0.25
bbox[4] <- max(w0$geometry[[1]][[1]][,2]) + 0.25

direction_df <- data.frame(xstart=double(), ystart=double(), xend=double(), yend=double(),
                           length=double(), angle=double())

for(i in 1:length(Lightning_data_list)){
  fit<-principal_curve(data.matrix(Lightning_data_list[[i]][,c("long", "lat", "tcount")]))

  curve_df <- as.data.frame(fit$s)
  colnames(curve_df) <- c("x", "y", "t")
  xstart=curve_df$x[1]
  ystart=curve_df$y[1]
  xend=curve_df$x[length(curve_df$x)]
  yend=curve_df$y[length(curve_df$y)] 
  arrow_length <- sqrt((xend-xstart)^2 + (yend-ystart)^2)
  arrow_rad <- atan2(yend-ystart, xend-xstart) 

  arrow_deg <- (arrow_rad * 180 / pi) %% 360
  arrow_deg
  direction_df_a <- data.frame(xstart=curve_df$x[1], ystart=curve_df$y[1],
                               xend=curve_df$x[length(curve_df$x)], yend=curve_df$y[length(curve_df$y)],
                               length=arrow_length, angle=arrow_deg)
  direction_df <- rbind(direction_df, direction_df_a)
                          
  #curve_df <- curve_df[order(fit$lambda), ]
  map <-get_stadiamap(bbox, maptype = "stamen_terrain_background", zoom=5) 
  p <- ggmap(map) + geom_sf(data=w0, inherit.aes = FALSE, fill=NA, color="black")  + geom_point(data=Lightning_data_list[[i]], aes(x=long,y=lat, color=tcount)) + geom_path(data = curve_df, aes(x = x, y = y), color="black") + scale_color_gradientn(colours = rainbow(20)) + labs(title=paste0("Group ", i, ", PC"), x = "Longitude", y = "Latitude")
  #ggsave(paste0("PC-",i,".pdf"), plot=p, width=5, height=3)
}

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


plot_gif <- FALSE
########################################
##Read data
########################################
Lightning_total <- readRDS("~/Library/CloudStorage/Dropbox/Data/Lightning/Lightning_2022.RDS")
#Lightning <- readRDS("~/Library/CloudStorage/Dropbox/Data/Lightning/Lightningevent.RDS")
Lightning <- readRDS("~/Library/CloudStorage/Dropbox/Data/Lightning/Lightningevent_2022.RDS")
#Lightning <- readRDS("C:\\Users\\HYU\\Dropbox\\Data\\Lightning\\Lightningevent.RDS")
unique_date <- unique(Lightning$data$일시)
#plot(Lightning$data$경도[which(Lightning$data$일시==unique_date[3])], Lightning$data$위도[which(Lightning$data$일시==unique_date[3])])


nrow(Lightning$data)/nrow(Lightning_total$data)

sum(month(Lightning_total$data$일시)%in%c(6,7,8))/nrow(Lightning_total$data)

table(cbind(month(Lightning_total$data$일시)))

table(cbind(month(Lightning$data$일시)))

w0 <- Lightning$chull

a <- Lightning_total$data[which(Lightning_total$data$일시>="2022-08-08 00:00:00 KST" & Lightning_total$data$일시<="2022-08-09 05:00:00 KST"),]
table(as.character(a$일시))

b <- Lightning$data[which(Lightning$data$일시>="2022-08-08 00:00:00 KST" & Lightning$data$일시<="2022-08-09 05:00:00 KST"),]
table(as.character(b$일시))

#plot(w0, axes = TRUE, col = 'white')
plot(w0)
maps::map(database = 'world',
          #regions = "antarctica",
          xlim = range(Lightning$data$경도),
          ylim = range(Lightning$data$위도),
          fill = T,
          col = 'white',
          resolution = 0,
          bg = 'white',
          mar = c(1,1,2,1), add=T
)
points(Lightning$data$경도, Lightning$data$위도, cex = .5, pch=16, col='red')

time_index_second <- seq(from = as.POSIXct("2019-01-01 00:00:00", tz="Asia/Seoul"), to = as.POSIXct("2023-12-31 23:59:59", tz="Asia/Seoul"), by = "sec")

#Lightning_time <-  floor_date(Lightning$data$시분초, "hour") 
Lightning$data$일시 <- as.POSIXct(Lightning$data$일시)
Lightning_count <- table(Lightning$data$일시)
plot(Lightning_count>50, type='l') #50회 이상

result_mat3 <- Lightning$info
result_mat3$nevent <- table(Lightning$data$event)
n_event <- length(unique(Lightning$data$event))

sum(result_mat3$duration)/(365*24) #총 5%정도
########################################
##plotting
########################################

order_size_vec <- c(11, 14, 2, 32, 33, 28, 12, 9, 31, 29, 10, 16, 22, 18, 26, 13, 34, 15, 8, 7, 19, 5, 17, 1, 4, 27, 3, 21,
                    24, 20, 25, 23, 30, 6)
library(ggplot2)
library(ggmap)
#NEED google keys
place <- "cheongju"
google <- get_googlemap(place, zoom = 10)
ggmap(google)

bbox <- bb2bbox(attr(google, "bb"))
bbox[1] <- min(w0$geometry[[1]][[1]][,1]) - 0.25
bbox[2] <- min(w0$geometry[[1]][[1]][,2]) - 0.25
bbox[3] <- max(w0$geometry[[1]][[1]][,1]) + 0.25
bbox[4] <- max(w0$geometry[[1]][[1]][,2]) + 0.25

map <-get_stadiamap(bbox, maptype = "stamen_terrain_background", zoom=5) 
ggmap(map) + geom_sf(data=w0, inherit.aes = FALSE, fill=NA, color="black") + labs(title="(a) Populations", x = "Longitude", y = "Latitude")

########################################
##summary statistics
########################################
#(1) year
table(year(Lightning$info$start))

#(2) month
table(month(Lightning$info$start))
hist(month(Lightning$info$start))

summary(Lightning$info$duration)
summary(as.numeric(table(Lightning$data$event)))

#pdf("duration.pdf", 10, 5)
#png("duration.png", 1000, 500)
par(mfrow=c(1,2))
hist(Lightning$info$duration, freq = T, xlab="Duration (hours)", main="(a) Histogram of Event Duration")
hist(table(Lightning$data$event), freq = T, xlab="No. of Lightnings", main="(b) Histogram of No. of Lightnings")
dev.off()

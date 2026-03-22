########################################
##Packages
########################################
library(spacetime)
library(gstat)
unloadNamespace(c('raster', 'terra')) 
require(xts)
require(forecast)
#for maps
require(maps)    # Provides functions that let us plot the maps
require(mapdata)    # Contains the hi-resolution points that mark out the countries.
require(sp)
require(rgdal)
require(maptools)
require(rgeos)
require(GeoXp)
require(ggplot2) #to use fortify function
library(shapefiles) #to simplyfy shape files
##For kriging
require(geoR) #for kriging
require(rrcov)
require(automap)
#for spatial GEVs
require(SpatialExtremes)
require(fExtremes) #for GEV estimation
require(extRemes) #for GEV estimation and model checking
require(ismev)
library(animation)
library(tidyverse)

library(lubridate)
library(sf)

library(PtProcess)
plot_gif <- FALSE
########################################
##Read data
########################################
#Lightning <- readRDS("~/Dropbox/Data/Lightning/Lightning_2019_2023.RDS")
Lightning <- readRDS("~/Dropbox/Data/Lightning/Lightning_2022.RDS")
#Lightning <- readRDS("C:\\Users\\HYU\\Dropbox\\Data\\Lightning\\Lightning_2019_2023.RDS")
unique_date <- unique(Lightning$data$일시)
#plot(Lightning$data$경도[which(Lightning$data$일시==unique_date[3])], Lightning$data$위도[which(Lightning$data$일시==unique_date[3])])

#plot(Lightning$data$경도, Lightning$data$위도)
#Observation window
#https://r-spatial.org/book/11-PointPattern.html
xx <- mean(c(min(Lightning$data$경도), max(Lightning$data$경도)))
yy <- mean(c(min(Lightning$data$위도), max(Lightning$data$위도)))
rr <- max(c(max(Lightning$data$경도)-xx, max(Lightning$data$위도)-yy))+1
ch <- chull(cbind(Lightning$data$경도, Lightning$data$위도))
ch_df <- data.frame(x=Lightning$data$경도[c(ch,ch[1])],
               y=Lightning$data$위도[c(ch,ch[1])])
ch0 <- ch_df %>%
  st_as_sf(coords=c('x','y')) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

w0 <- st_sfc(st_point(c(xx, yy-1))) |> st_buffer(rr)
plot(w0, axes = TRUE, col = 'grey')
points(Lightning$data$경도, Lightning$data$위도, add = TRUE, cex = .5, pch=16, col='red')


time_index_second <- seq(from = as.POSIXct("2022-01-01 00:00:00", tz="Asia/Seoul"), to = as.POSIXct("2022-12-31 23:59:59", tz="Asia/Seoul"), by = "sec")

#Lightning_time <-  floor_date(Lightning$data$시분초, "hour") 
Lightning$data$일시 <- as.POSIXct(Lightning$data$일시)
Lightning_count <- table(Lightning$data$일시)
plot(Lightning_count>50, type='l') 

#data <- Lightning$data; threshold <- 50
event_detection <- function(data, threshold=50, duration=NA) {
  Lightning_count <- table(data$일시)
  Lightning_count_name <- names(Lightning_count)
  Lightning_POS <- as.POSIXct(Lightning_count_name, tz="Asia/Seoul", format=c("%Y-%m-%d %H:%M:%OS"))
  Lightning_POS[which(is.na(Lightning_POS))] <- format(as.POSIXct(Lightning_count_name, tz="Asia/Seoul", format=c("%Y-%m-%d")),"%Y-%m-%d %H:%M:%OS")[which(is.na(Lightning_POS))]
  
  which_index <- which(Lightning_count>threshold)
  Lightning_count <- Lightning_count[which_index]
  Lightning_POS <- Lightning_POS[which_index]
  Lightning_POS[which(is.na(Lightning_POS))] <- format(as.POSIXct(Lightning_count_name, tz="Asia/Seoul", format=c("%Y-%m-%d")),"%Y-%m-%d %H:%M:%OS")[which(is.na(Lightning_POS))]
  result_mat <- c()
  duration_vec <- c()
  n <- length(Lightning_POS)
  i <- 1
  while(TRUE){
    #cat(i)
    newind <- c(Lightning_POS[i])
    duration <- 0
    while(TRUE){
      if(i==n){
        newind <- c(newind, Lightning_POS[i])
        break
      }else{
        compind <- Lightning_POS[i] + lubridate::minutes(60)
        if(compind==Lightning_POS[(i+1)]){
          i <- (i+1)
          duration <- duration + 1
        }else{
          newind <- c(newind, Lightning_POS[i])
          #print(newind)
          i <- (i+1)
          duration <- duration + 1
          break
        }
      }
    }
    result_mat <- rbind(result_mat, newind)
    duration_vec <- c(duration_vec, duration)
    if(i >= n){
      break
    }
  }
  return(data.frame(start=as.POSIXct(result_mat[,1]), end=as.POSIXct(result_mat[,2])+ lubridate::minutes(60), duration=duration_vec))
}

result_mat2 <- event_detection(data = Lightning$data, threshold = 50)
#즉 start, (원래) end time이 2/14 19:00으로 같으면
#2/14 19:00 ~ 2/14 19:59 낙뢰 이벤트가 50개 이상이라는 뜻이므로 duration이 1시간
#end time에 1시간 추가해서 저장
plot(result_mat2$duration)
hist(result_mat2$duration)
quantile(result_mat2$duration, probs=0.75) #8
quantile(result_mat2$duration, probs=0.8) #9
quantile(result_mat2$duration, probs=0.85) #12
quantile(result_mat2$duration, probs=0.9) #16
sum(result_mat2$duration>=8) #202

result_mat3 <- result_mat2[which(result_mat2$duration>=8),]

if(plot_gif==TRUE){
  par(mfrow=c(1,1))
  #plot_index <- 3
  for(plot_index in 1:length(result_mat3$duration)){
    col_palette <- rainbow(result_mat3$duration[plot_index])
    saveGIF({
      for(i in 1:result_mat3$duration[plot_index]){
        current_time <- result_mat3$start[plot_index] + lubridate::minutes(60*(i-1))
        #Lightning$data$일시
        plot(Lightning$shape, main=current_time
             , xlim=range(Lightning$data$경도), ylim=range(Lightning$data$위도))
        points(Lightning$data$경도[which(Lightning$data$일시==current_time)], Lightning$data$위도[which(Lightning$data$일시==current_time)], pch=16, col=col_palette[i])
      }
    }, movie.name = paste0("Lightning",plot_index,".gif"),
    ani.width=700, ani.height=700)
  }
}


########################################
##Save 34 events
########################################
Lightning.data2 <- data.frame()
for(i in 1:nrow(result_mat3)){
  interval1 <- interval(result_mat3$start[i]-hours(1), result_mat3$end[i]+hours(1))
  Lightning.data2 <- rbind(Lightning.data2, cbind(Lightning$data[which(Lightning$data$시분초 %within% interval1),], event=as.factor(i)))
}

dim(Lightning.data2)/dim(Lightning$data)

xx <- mean(c(min(Lightning$data$경도), max(Lightning$data$경도)))
yy <- mean(c(min(Lightning$data$위도), max(Lightning$data$위도)))
rr <- max(c(max(Lightning$data$경도)-xx, max(Lightning$data$위도)-yy))+1

w0 <- st_sfc(st_point(c(xx, yy-1))) |> st_buffer(rr)
plot(w0, axes = TRUE, col = 'grey')


if(plot_gif==TRUE){
  par(mfrow=c(1,1))
  #plot_index <- 3
  for(plot_index in c(10,16)){
    col_palette <- rainbow(2)
    saveGIF({
      for(i in 1:(6*result_mat3$duration[plot_index]-5)){
        current_time0 <- result_mat3$start[plot_index] + lubridate::minutes(10*(i-1)-10)
        current_time <- result_mat3$start[plot_index] + lubridate::minutes(10*(i-1))
        current_time2 <- result_mat3$start[plot_index] + lubridate::minutes(10*(i-1)+10)
        #current_time3 <- result_mat3$start[plot_index] + lubridate::minutes(10*(i-1)+20)
        #current_time4 <- result_mat3$start[plot_index] + lubridate::minutes(60*(i-1)+30)
        #current_time5 <- result_mat3$start[plot_index] + lubridate::minutes(60*(i-1)+40)
        #current_time6 <- result_mat3$start[plot_index] + lubridate::minutes(60*(i-1)+50)
        #Lightning$data$일시
        plot(Lightning$shape, main=current_time
             , xlim=range(Lightning$data$경도[which(Lightning.data2$event==plot_index)]), ylim=range(Lightning$data$위도[which(Lightning.data2$event==plot_index)]))
        points(Lightning$data$경도[which(Lightning$data$시분초>=current_time0  & Lightning$data$시분초<=current_time)], Lightning$data$위도[which(Lightning$data$시분초>=current_time0  & Lightning$data$시분초<=current_time)], pch=16, col=col_palette[1])
        points(Lightning$data$경도[which(Lightning$data$시분초>=current_time  & Lightning$data$시분초<=current_time2)], Lightning$data$위도[which(Lightning$data$시분초>=current_time  & Lightning$data$시분초<=current_time2)], pch=16, col=col_palette[2])
        #points(Lightning$data$경도[which(Lightning$data$시분초>=current_time2  & Lightning$data$시분초<=current_time3)], Lightning$data$위도[which(Lightning$data$시분초>=current_time2  & Lightning$data$시분초<=current_time3)], pch=16, col=col_palette[3])
        #points(Lightning$data$경도[which(Lightning$data$시분초>=current_time3  & Lightning$data$시분초<=current_time4)], Lightning$data$위도[which(Lightning$data$시분초>=current_time3  & Lightning$data$시분초<=current_time4)], pch=16, col=col_palette[4])
        #points(Lightning$data$경도[which(Lightning$data$시분초>=current_time4  & Lightning$data$시분초<=current_time5)], Lightning$data$위도[which(Lightning$data$시분초>=current_time4  & Lightning$data$시분초<=current_time5)], pch=16, col=col_palette[5])
        #points(Lightning$data$경도[which(Lightning$data$시분초>=current_time5  & Lightning$data$시분초<=current_time6)], Lightning$data$위도[which(Lightning$data$시분초>=current_time5  & Lightning$data$시분초<=current_time6)], pch=16, col=col_palette[6])
      }
    }, movie.name = paste0("Lightning",plot_index,".gif"),
    ani.width=700, ani.height=700)
  }
}

Lightning_list <- list()
Lightning_list$data <- Lightning.data2
Lightning_list$shape <- Lightning$shape
Lightning_list$info <- result_mat3
Lightning_list$window <- w0
Lightning_list$chull <- ch0
saveRDS(Lightning_list, "~/Library/CloudStorage/Dropbox/Data/Lightning/Lightningevent_2022.RDS")


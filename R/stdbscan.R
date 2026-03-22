########################################
##Packages
########################################
library(stdbscan)
library(readr)
library(ggplot2)
library(plotly)
library(xts)
library(lubridate)
########################################
##Read data
########################################
Lightning <- readRDS("~/Dropbox/Data/Lightning/Lightning_2022.RDS")
#Lightning <- readRDS("~/Library/CloudStorage/Dropbox/Data/Lightning/Lightning_2022.RDS")


########################################
##Preprocessing
##https://cran.r-project.org/web/packages/stdbscan/vignettes/stop-identification.html
##ST-DBSCAN https://www.sciencedirect.com/science/article/pii/S0169023X06000218
########################################
##(EPSG:4586)
library(sf)
# (WGS84)
df_sf <- st_as_sf(Lightning$data, coords = c("경도", "위도"), crs = 4326)
# EPSG:4586
df_sf_proj <- st_transform(df_sf, 4586)
coords <- st_coordinates(df_sf_proj)

t0 <- as.POSIXct("2022-01-01 00:00:00", tz = "Asia/Seoul")
sec <- as.numeric(difftime(index(Lightning$data), t0, units = "secs"))

(res <- st_dbscan(
  x=coords[,1], y=coords[,2], t=sec,
  eps_spatial = 1000*20, # meters
  eps_temporal = 60*60*8, # seconds
  min_pts = 50*8
))

max(res)

library(dplyr)

res_list <- list()

#manual search

res1 <- st_dbscan(
  x = coords[,1],
  y = coords[,2],
  t = sec,
  eps_spatial  = 16000,
  eps_temporal = 28800 , 
  min_pts      = 400
)

plot(res1)


cluster_ids <- c(1:34)

result_stdbscan1 <- data.frame()
for (cid in cluster_ids) {
  times <- Lightning$data$일시[which(res1==cid)]
  
  start_time <- min(times)
  end_time   <- max(times)
  
  duration <- as.numeric(difftime(end_time, start_time, units = "hours"))
  
  result_stdbscan1 <- rbind(result_stdbscan1, data.frame(
    #cluster_id = cid,
    start_time = start_time,
    end_time   = end_time,
    duration_hours = duration
  ))
}



res2 <- st_dbscan(
  x = coords[,1],
  y = coords[,2],
  t = sec,
  eps_spatial  = 16000*2,
  eps_temporal = 2880/10 , 
  min_pts      = 410
)
plot(res2)
plot(res2[which(res2!=-1)])
max(res2)

result_stdbscan2 <- data.frame()
for (cid in cluster_ids) {
  times <- Lightning$data$일시[which(res2==cid)]
  
  start_time <- min(times)
  end_time   <- max(times)
  
  duration <- as.numeric(difftime(end_time, start_time, units = "hours"))
  
  result_stdbscan2 <- rbind(result_stdbscan2, data.frame(
    #cluster_id = cid,
    start_time = start_time,
    end_time   = end_time,
    duration_hours = duration
  ))
}


res3 <- st_dbscan(
  x = coords[,1],
  y = coords[,2],
  t = sec,
  eps_spatial  = 16000*5.4,
  eps_temporal = 2880/1.75 , 
  min_pts      = 750
)
plot(res3)
plot(res3[which(res2!=-1)])
max(res3)


length(Lightning$data$일시[which(res2!=-1)])/length(Lightning$data$일시)

cluster_ids <- c(1:34)

result_stdbscan3 <- data.frame()
for (cid in cluster_ids) {
  times <- Lightning$data$일시[which(res3==cid)]
  
  start_time <- min(times)
  end_time   <- max(times)
  
  duration <- as.numeric(difftime(end_time, start_time, units = "hours"))
  
  result_stdbscan3 <- rbind(result_stdbscan3, data.frame(
    #cluster_id = cid,
    start_time = start_time,
    end_time   = end_time,
    duration_hours = duration
  ))
}
result_stdbscan3

res_list[[1]] <- res1
res_list[[2]] <- res2
res_list[[3]] <- res3
res_list[[4]] <- result_stdbscan1
res_list[[5]] <- result_stdbscan2
res_list[[6]] <- result_stdbscan3
res_list[[7]] <- data.frame(names=c("res1","res2","res3"),
                            eps_spatial  = c(16000, 16000*2, 16000*5.4),
                            eps_temporal = c(28800, 2880/10, 2880/1.75),
                            min_pts      = c(400,410,750))

saveRDS(res_list, "rev_stdbscan.RDS")

#res_list <- readRDS("~/Dropbox/R files/Lightning/revision/rev_stdbscan.RDS")
#res1 <- res_list[[1]]; res2 <- res_list[[2]]; res3 <- res_list[[3]]; result_stdbscan1 <- res_list[[4]]; result_stdbscan2 <- res_list[[5]]; result_stdbscan3 <- res_list[[6]]
#Lightningevent <- readRDS("~/Dropbox/Data/Lightning/Lightningevent_2022.RDS")
#Lightningevent <- readRDS("~/Library/CloudStorage/Dropbox/Data/Lightning/Lightningevent_2022.RDS")
Lightningevent$info #for comparison


#cluster coverage
nrow(Lightningevent$data)/nrow(Lightning$data)
length(Lightning$data$일시[which(res1!=-1)])/length(Lightning$data$일시)
length(Lightning$data$일시[which(res2!=-1)])/length(Lightning$data$일시)
length(Lightning$data$일시[which(res3!=-1)])/length(Lightning$data$일시)

#median of cluster size
median(table(Lightningevent$data$event))
median(table(res1)[-1])
median(table(res2)[-1])
median(table(res3)[-1])

#median of cluster duration
median(Lightningevent$info$duration)
median(result_stdbscan1$duration_hours)
median(result_stdbscan2$duration_hours)
median(result_stdbscan3$duration_hours)

#mean of cluster duration
mean(Lightningevent$info$duration)
mean(result_stdbscan1$duration_hours)
mean(result_stdbscan2$duration_hours)
mean(result_stdbscan3$duration_hours)

#maximum of cluster duration
max(Lightningevent$info$duration)
max(result_stdbscan1$duration_hours)
max(result_stdbscan2$duration_hours)
max(result_stdbscan3$duration_hours)

#Jaccard index
key_event <- paste(
  Lightningevent$data$시분초,
  Lightningevent$data$경도,
  Lightningevent$data$위도
)

key_all <- paste(
  Lightning$data$시분초,
  Lightning$data$경도,
  Lightning$data$위도
)

key_res1 <- paste(
  Lightning$data$시분초[which(res1!=-1)],
  Lightning$data$경도[which(res1!=-1)],
  Lightning$data$위도[which(res1!=-1)]
)

key_res2 <- paste(
  Lightning$data$시분초[which(res2!=-1)],
  Lightning$data$경도[which(res2!=-1)],
  Lightning$data$위도[which(res2!=-1)]
)

key_res3 <- paste(
  Lightning$data$시분초[which(res3!=-1)],
  Lightning$data$경도[which(res3!=-1)],
  Lightning$data$위도[which(res3!=-1)]
)

jaccard_key <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  length(intersect(a, b)) / length(union(a, b))
}


j_event_res1 <- jaccard_key(key_event, key_res1)
j_event_res2 <- jaccard_key(key_event, key_res2)
j_event_res3 <- jaccard_key(key_event, key_res3)

j_event_res1
j_event_res2
j_event_res3


pr_metrics <- function(truth, pred) {
  truth <- unique(truth)
  pred  <- unique(pred)
  
  TP <- length(intersect(truth, pred))
  FP <- length(setdiff(pred, truth))
  FN <- length(setdiff(truth, pred))
  
  precision <- TP / (TP + FP)
  recall    <- TP / (TP + FN)
  f1        <- 2 * TP / (2 * TP + FP + FN)
  jaccard   <- TP / (TP + FP + FN)
  
  c(
    precision = precision,
    recall = recall,
    f1 = f1,
    jaccard = jaccard
  )
}

res1_metrics <- pr_metrics(key_event, key_res1)
res2_metrics <- pr_metrics(key_event, key_res2)
res3_metrics <- pr_metrics(key_event, key_res3)

metrics_df <- rbind(
  "ST-DBSCAN 1" = res1_metrics,
  "ST-DBSCAN 2" = res2_metrics,
  "ST-DBSCAN 3" = res3_metrics
)

metrics_df
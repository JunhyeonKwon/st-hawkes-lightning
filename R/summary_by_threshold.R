########################################
##Goal
##(1) number of clusters
##(2) median size
##(3) median duration (hours)
########################################
#/CloudStorage/Dropbox
Lightning_2022_40 <- readRDS("~/Dropbox/Data/Lightning/Lightningevent_2022_40.RDS")
Lightning_2022_45 <- readRDS("~/Dropbox/Data/Lightning/Lightningevent_2022_45.RDS")
Lightning_2022_50 <- readRDS("~/Dropbox/Data/Lightning/Lightningevent_2022.RDS")
Lightning_2022_55 <- readRDS("~/Dropbox/Data/Lightning/Lightningevent_2022_55.RDS")
Lightning_2022_60 <- readRDS("~/Dropbox/Data/Lightning/Lightningevent_2022_60.RDS")

Lightning_2023_40 <- readRDS("~/Dropbox/Data/Lightning/Lightningevent_2023_40.RDS")
Lightning_2023_45 <- readRDS("~/Dropbox/Data/Lightning/Lightningevent_2023_45.RDS")
Lightning_2023_50 <- readRDS("~/Dropbox/Data/Lightning/Lightningevent_2023.RDS")
Lightning_2023_55 <- readRDS("~/Dropbox/Data/Lightning/Lightningevent_2023_55.RDS")
Lightning_2023_60 <- readRDS("~/Dropbox/Data/Lightning/Lightningevent_2023_60.RDS")
########################################
##(1) number of clusters
########################################
length(unique(Lightning_2022_40$data$event))
length(unique(Lightning_2022_45$data$event))
length(unique(Lightning_2022_50$data$event))
length(unique(Lightning_2022_55$data$event))
length(unique(Lightning_2022_60$data$event))

length(unique(Lightning_2023_40$data$event))
length(unique(Lightning_2023_45$data$event))
length(unique(Lightning_2023_50$data$event))
length(unique(Lightning_2023_55$data$event))
length(unique(Lightning_2023_60$data$event))

########################################
##(2) median size
########################################
median(table(Lightning_2022_40$data$event))
median(table(Lightning_2022_45$data$event))
median(table(Lightning_2022_50$data$event))
median(table(Lightning_2022_55$data$event))
median(table(Lightning_2022_60$data$event))

median(table(Lightning_2023_40$data$event))
median(table(Lightning_2023_45$data$event))
median(table(Lightning_2023_50$data$event))
median(table(Lightning_2023_55$data$event))
median(table(Lightning_2023_60$data$event))

########################################
##(3) median duration (hours)
########################################
median(Lightning_2022_40$info$duration)
median(Lightning_2022_45$info$duration)
median(Lightning_2022_50$info$duration)
median(Lightning_2022_55$info$duration)
median(Lightning_2022_60$info$duration)

median(Lightning_2023_40$info$duration)
median(Lightning_2023_45$info$duration)
median(Lightning_2023_50$info$duration)
median(Lightning_2023_55$info$duration)
median(Lightning_2023_60$info$duration)
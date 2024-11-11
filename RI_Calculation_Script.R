####### Script to calculate RI #######

# Installation of R package to calculate linear retention index (RI).

# Installation of "MetaboCoreUtils" package
#install.packages("remotes")
#remotes::install_github("rformassspectrometry/MetaboCoreUtils")

# Loading "MetaboCoreUtils" library
library("MetaboCoreUtils")

# RI for in-house library

## Loadding rt of each n-alkane
### chromatogram 1 "2_n-Alkanes"
rtime_lib1 <- c(3.475, 5.105, 7.24, 9.685, 12.235, 14.76, 17.215, 19.56, 21.8,
                23.935, 25.975, 28.23, 31.185, 35.21, 39.225, 41.925, 44.04,
                45.825, 47.39, 48.81, 50.11, 51.335, 52.485, 53.645, 54.965)
rindex_lib1 <- c(900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
                 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800,
                 2900, 3000, 3100, 3200, 3300)
rti_lib1 <- data.frame(rtime = rtime_lib1, rindex = rindex_lib1)
### chromatogram 2 "7_n-Alkanes"
rtime_lib2 <- c(3.466, 5.095, 7.225, 9.67, 12.215, 14.745, 17.201, 19.546,
                21.785, 23.92, 25.96, 28.21, 31.165, 35.17, 39.195, 41.905,
                44.02, 45.811, 47.375, 48.796, 50.095, 51.325, 52.475, 53.635,
                54.945)
rti_lib2 <- data.frame(rtime = rtime_lib2, rindex = rindex_lib1)
### chromatogram 3 "14_n-Alkanes"
rtime_lib3 <- c(3.47, 5.101, 7.235, 9.68, 12.231, 14.755, 17.215, 19.555,
                21.796, 23.93, 25.971, 28.226, 31.18, 35.195, 39.215, 41.92,
                44.04, 45.82, 47.391, 48.806, 50.11, 51.335, 52.485, 53.645,
                54.961)
rti_lib3 <- data.frame(rtime = rtime_lib3, rindex = rindex_lib1)
### chromatogram 4 "19_n-Alkanes"
rtime_lib4 <- c(3.466, 5.096, 7.23, 9.675, 12.22, 14.75, 17.21, 19.545, 21.79,
                23.925, 25.965, 28.221, 31.17, 35.19, 39.21, 41.915, 44.03,
                45.816, 47.385, 48.8, 50.1, 51.325, 52.475, 53.635, 54.946)
rti_lib4 <- data.frame(rtime = rtime_lib4, rindex = rindex_lib1)
### chromatogram 5 "26_n-Alkanes"
rtime_lib5 <- c(3.465, 5.095, 7.23, 9.675, 12.22, 14.751, 17.205, 19.546,
                21.785, 23.926, 25.96, 28.21, 31.165, 35.17, 39.2, 41.905,
                44.025, 45.81, 47.38, 48.795, 50.095, 51.32, 52.475, 53.63,
                54.945)
rti_lib5 <- data.frame(rtime = rtime_lib5, rindex = rindex_lib1)
### chromatogram 6 "31_n-Alkanes"
rtime_lib6 <- c(3.47, 5.095, 7.23, 9.675, 12.22, 14.75, 17.21, 19.545, 21.79,
                23.925, 25.965, 28.22, 31.175, 35.19, 39.215, 41.911, 44.035,
                45.815, 47.38, 48.8, 50.1, 51.33, 52.48, 53.64, 54.95)
rti_lib6 <- data.frame(rtime = rtime_lib6, rindex = rindex_lib1)
## Experimental Retention Index (RI) calculation
##with chromatogram 1 "2_n-Alkanes"
indexRtime(rtime_lib2, rti_lib1) #chromatogram 2
indexRtime(rtime_lib3, rti_lib1) #chromatogram 3
indexRtime(rtime_lib4, rti_lib1) #chromatogram 4
indexRtime(rtime_lib5, rti_lib1) #chromatogram 5
indexRtime(rtime_lib6, rti_lib1) #chromatogram 6
##with chromatogram 2 "7_n-Alkanes"
indexRtime(rtime_lib1, rti_lib2) #chromatogram 1
indexRtime(rtime_lib3, rti_lib2) #chromatogram 3
indexRtime(rtime_lib4, rti_lib2) #chromatogram 4
indexRtime(rtime_lib5, rti_lib2) #chromatogram 5
indexRtime(rtime_lib6, rti_lib2) #chromatogram 6
##with chromatogram 3 "14_n-Alkanes"
indexRtime(rtime_lib1, rti_lib3) #chromatogram 1
indexRtime(rtime_lib2, rti_lib3) #chromatogram 2
indexRtime(rtime_lib4, rti_lib3) #chromatogram 4
indexRtime(rtime_lib5, rti_lib3) #chromatogram 5
indexRtime(rtime_lib6, rti_lib3) #chromatogram 6
##with chromatogram 4 "19_n-Alkanes"
indexRtime(rtime_lib1, rti_lib4) #chromatogram 1
indexRtime(rtime_lib2, rti_lib4) #chromatogram 2
indexRtime(rtime_lib3, rti_lib4) #chromatogram 3
indexRtime(rtime_lib5, rti_lib4) #chromatogram 5
indexRtime(rtime_lib6, rti_lib4) #chromatogram 6
##with chromatogram 5 "26_n-Alkanes"
indexRtime(rtime_lib1, rti_lib5) #chromatogram 1
indexRtime(rtime_lib2, rti_lib5) #chromatogram 2
indexRtime(rtime_lib3, rti_lib5) #chromatogram 3
indexRtime(rtime_lib4, rti_lib5) #chromatogram 4
indexRtime(rtime_lib6, rti_lib5) #chromatogram 6
##with chromatogram 6 "31_n-Alkanes"
indexRtime(rtime_lib1, rti_lib6) #chromatogram 1
indexRtime(rtime_lib2, rti_lib6) #chromatogram 2
indexRtime(rtime_lib3, rti_lib6) #chromatogram 3
indexRtime(rtime_lib4, rti_lib6) #chromatogram 4
indexRtime(rtime_lib5, rti_lib6) #chromatogram 5








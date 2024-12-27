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

# Experimental RI for features that match with NIST libraries

## Loadding rt of each n-alkane
rtime_ids <- c(7.557, 10.006, 12.569, 15.111, 17.581, 19.937, 22.190, 24.338,
               26.399, 28.813, 32.000, 36.327, 39.949, 42.506, 44.557, 46.309,
               47.852, 49.257, 50.554, 51.781, 52.936, 54.182, 55.604)
rindex_ids <- c(1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
                2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
                3100, 3200, 3300)
rti_ids <- data.frame(rtime = rtime_ids, rindex = rindex_ids)
## Experimental Retention Index (RI) calculation
##Feature deconvoluted with MZmine
rtime_mzmine <- c(4.107, 4.293, 5.341, 6.077, 6.098, 6.772, 8.002, 8.376, 9.121,
                  9.971, 10.474, 10.835, 11.571, 11.845, 11.860, 12.310, 12.490,
                  12.609, 13.004, 13.205, 13.472, 13.924, 14.296, 14.706,
                  14.938, 16.919, 17.811, 21.068, 21.152, 22.327, 22.928,
                  23.786, 25.633, 25.698, 25.832, 25.835, 25.849, 25.851,
                  26.409, 26.526, 26.880, 26.956, 27.801, 27.848, 27.852,
                  27.953, 28.551, 28.634, 32.468, 33.543, 39.918, 42.484,
                  44.541, 49.664)
indexRtime(rtime_mzmine, rti_ids)
##Feature deconvoluted with eRah
rtime_erah <- c(4.079, 4.081, 4.082, 4.086, 4.088, 4.292, 4.295, 4.297, 4.297,
                5.130, 5.131, 5.294, 5.301, 5.306, 6.442, 6.442, 7.433, 7.481,
                7.484, 7.982, 7.987, 7.988, 7.989, 8.000, 8.265, 8.368, 9.050,
                9.066, 9.765, 9.767, 9.967, 9.969, 10.466, 10.473, 10.824,
                10.967, 11.140, 11.819, 12.485, 12.489, 12.496, 12.594, 12.999,
                13.001, 13.200, 13.203, 13.203, 13.449, 13.449, 13.450, 13.640,
                13.920, 13.920, 13.922, 13.923, 14.289, 14.759, 14.762, 15.702,
                16.306, 16.306, 16.906, 16.909, 16.915, 16.923, 19.173, 19.548,
                19.551, 20.663, 21.069, 21.069, 21.156, 21.337, 21.337, 21.504,
                21.667, 21.788, 21.974, 22.268, 22.268, 22.327, 22.329, 22.733,
                22.733, 22.770, 22.773, 22.860, 22.861, 22.979, 23.306, 23.768,
                23.786, 23.788, 24.522, 24.523, 25.424, 25.426, 25.537, 25.669,
                25.711, 25.850, 25.976, 26.035, 26.290, 26.290, 26.291, 26.519,
                26.604, 26.961, 27.253, 27.257, 27.742, 27.745, 27.847, 27.850,
                27.851, 28.557, 28.638, 28.640, 31.996, 32.018, 39.100, 39.102,
                39.104, 39.923, 39.923, 42.486, 42.487, 44.542, 44.543, 47.839,
                48.954, 49.668)
indexRtime(rtime_erah, rti_ids)
##Feature deconvoluted with MSHub
rtime_mshub <- c(3.140, 3.590, 3.750, 3.880, 3.950, 4.100, 4.250, 4.250, 4.310,
                 4.690, 4.770, 4.770, 4.920, 4.920, 5.000, 5.160, 5.300, 5.510,
                 5.570, 5.820, 6.220, 6.290, 6.400, 6.450, 6.650, 6.810, 6.950,
                 7.330, 7.430, 7.430, 7.600, 7.710, 7.880, 7.880, 8.000, 8.090,
                 8.160, 8.270, 8.270, 8.370, 8.610, 8.710, 9.020, 9.020, 9.160,
                 9.160, 9.340, 9.440, 9.440, 9.520, 9.670, 9.670, 9.830, 9.970,
                 10.280, 10.370, 10.470, 10.950, 11.060, 11.570, 11.830, 12.070,
                 12.330, 12.480, 12.600, 12.750, 12.900, 13.000, 13.000, 13.200,
                 13.440, 13.650, 13.650, 13.840, 13.910, 14.290, 14.550, 14.690,
                 15.140, 15.450, 15.700, 16.580, 16.720, 16.840, 17.660, 17.870,
                 18.150, 18.350, 19.020, 19.170, 19.380, 19.540, 20.340, 20.590,
                 20.660, 21.330, 21.980, 22.260, 22.320, 22.550, 22.730, 22.850,
                 22.900, 22.970, 23.230, 23.780, 24.320, 24.520, 25.210, 25.630,
                 25.740, 25.840, 25.840, 25.970, 26.030, 26.280, 26.540, 26.950,
                 27.120, 27.480, 27.740, 27.850, 27.850, 27.950, 28.540, 28.630,
                 28.630, 28.790, 29.140, 30.910, 31.770, 32.000, 32.000, 32.210,
                 32.210, 32.460, 33.200, 3.530, 34.010, 34.320, 34.460, 34.770,
                 34.950, 36.030, 36.120, 36.290, 36.750, 37.200, 37.460, 37.510,
                 37.830, 38.170, 38.490, 38.890, 39.090, 39.310, 39.920, 41.370,
                 42.380, 42.480, 42.700, 43.100, 44.010, 44.210, 44.270, 44.370,
                 44.540, 44.650, 46.290, 47.090, 47.830, 48.950, 49.540, 49.660,
                 49.810, 50.530, 51.020, 51.750)
indexRtime(rtime_mshub, rti_ids)


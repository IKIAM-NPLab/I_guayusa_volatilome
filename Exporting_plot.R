
# Exporting plot

## General PCA

#Merge of score and loading plots of PCA analysis.



# Library loading
library(gridExtra)
# Figure matrix
figure_1 <- arrangeGrob(figure_1a,
                        figure_1b,
                        figure_1c,
                        figure_1d,
                        layout_matrix = rbind(c(1, 2, 3),
                                              c(4, 4, 4),
                                              c(4, 4, 4)))
# Adding label to the figures
figure_one <- ggpubr::as_ggplot(figure_1) +
  draw_plot_label(label = LETTERS[1:4],
                  x = c(0, 0.332, 0.665, 0),
                  y = c(.99, .99, .99, .656))
# Exporting (*.pdf) file
#ggsave(filename = "Result/notame_Result/figure_1_glog.pdf", plot = figure_one,
#      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
#ggsave(filename = "Result/notame_Result/figure_1_glog.png", plot = figure_one,
#      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.jpg) file
#ggsave(filename = "Result/notame_Result/figure_1_glog.jpg", plot = figure_one,
#      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)



#PCA score plots of different factors.



# Figure matrix
figure_1_pca <- arrangeGrob(figure_1a,
                            figure_1b,
                            figure_1c,
                            layout_matrix = rbind(c(1, 2, 3)))
# Adding label to the figures
figure_one_pca <- ggpubr::as_ggplot(figure_1_pca) +
  draw_plot_label(label = LETTERS[1:3],
                  x = c(0, 0.332, 0.665),
                  y = c(.99, .99, .99))
# Exporting (*.pdf) file
#ggsave(filename = "Result/notame_Result/figure_1_glog_pca.pdf",
#       plot = figure_one_pca, width = 175, height = 45, units = "mm",
#       dpi = 300, scale = 2.5)
# Exporting (*.png) file
#ggsave(filename = "Result/notame_Result/figure_1_glog_pca.png",
#       plot = figure_one_pca, width = 175, height = 45, units = "mm",
#       dpi = 300, scale = 2.5)
# Exporting (*.jpg) file
#ggsave(filename = "Result/notame_Result/figure_1_glog_pca.jpg",
#       plot = figure_one_pca, width = 175, height = 45, units = "mm",
#       dpi = 300, scale = 2.5)



## PCA by location

### Alto Pano



# Figure matrix
figure_s1 <- arrangeGrob(figure_s1a,
                         figure_s1b,
                         figure_s1c,
                         figure_s1d,
                         layout_matrix = rbind(c(1, 2, 3),
                                               c(4, 4, 4),
                                               c(4, 4, 4)))
# Adding label to the figures
figure_sone <- ggpubr::as_ggplot(figure_s1) +
  draw_plot_label(label = LETTERS[1:4],
                  x = c(0, 0.332, 0.665, 0),
                  y = c(.99, .99, .99, .656))
# Exporting (*.pdf) file
#ggsave(filename = "Result/notame_Result/figure_s1.pdf", plot = figure_sone,
#      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
#ggsave(filename = "Result/notame_Result/figure_s1.png", plot = figure_sone,
#      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.jpg) file
#ggsave(filename = "Result/notame_Result/figure_s1.jpg", plot = figure_sone,
#      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)



### Alto Tena



# Figure matrix
figure_s2 <- arrangeGrob(figure_s1a,
                         figure_s1b,
                         figure_s1c,
                         figure_s1e,
                         layout_matrix = rbind(c(1, 2, 3),
                                               c(4, 4, 4),
                                               c(4, 4, 4)))
# Adding label to the figures
figure_stwo <- ggpubr::as_ggplot(figure_s2) +
  draw_plot_label(label = LETTERS[1:4],
                  x = c(0, 0.332, 0.665, 0),
                  y = c(.99, .99, .99, .656))
# Exporting (*.pdf) file
#ggsave(filename = "Result/notame_Result/figure_s2.pdf", plot = figure_stwo,
#      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
#ggsave(filename = "Result/notame_Result/figure_s2.png", plot = figure_stwo,
#      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.jpg) file
#ggsave(filename = "Result/notame_Result/figure_s2.jpg", plot = figure_stwo,
#      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)



### Talag



# Figure matrix
figure_s3 <- arrangeGrob(figure_s1a,
                         figure_s1b,
                         figure_s1c,
                         figure_s1f,
                         layout_matrix = rbind(c(1, 2, 3),
                                               c(4, 4, 4),
                                               c(4, 4, 4)))
# Adding label to the figures
figure_sthree <- ggpubr::as_ggplot(figure_s3) +
  draw_plot_label(label = LETTERS[1:4],
                  x = c(0, 0.332, 0.665, 0),
                  y = c(.99, .99, .99, .656))
# Exporting (*.pdf) file
#ggsave(filename = "Result/notame_Result/figure_s3.pdf", plot = figure_sthree,
#      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
#ggsave(filename = "Result/notame_Result/figure_s3.png", plot = figure_sthree,
#      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.jpg) file
#ggsave(filename = "Result/notame_Result/figure_s3.jpg", plot = figure_sthree,
#      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)



## Age plots



# Figure matrix
figure_2 <- arrangeGrob(gcms_hm,
                        vc_plot,
                        layout_matrix = rbind(c(1, 2)))
# Adding label to the figures
figure_two <- ggpubr::as_ggplot(figure_2) +
  draw_plot_label(label = LETTERS[1:2],
                  x = c(0, 0.5),
                  y = c(.99, .99))
# Exporting (*.pdf) file
#ggsave(filename = "Result/notame_Result/figure_2.pdf", plot = figure_two,
#      width = 175, height = 155, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
#ggsave(filename = "Result/notame_Result/figure_2.png", plot = figure_two,
#      width = 175, height = 155, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.jpg) file
#ggsave(filename = "Result/notame_Result/figure_2.jpg", plot = figure_two,
#      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)



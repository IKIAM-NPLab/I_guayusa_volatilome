
# Heatmap for Journal of Food Composition and Analysis (JFCA)

# Extracting metabolites with statistic significance
hm_jfca <- hm_age_set[hm_age_set@featureData@data$Early != "ns",]
# Heatmap
set.seed(2571)
# Extracting feature height table
hm_height <- exprs(hm_jfca)
# Extracting sample information
hm_pdata <- hm_jfca@phenoData@data
# Extracting feature information
hm_fdata <- hm_jfca@featureData@data
# Adding row and column name
rownames(hm_height) <- hm_fdata$Metabolite
colnames(hm_height) <- hm_pdata$Group
# Metabolite superclass color
cols_metclass <- c("Organic oxygen compounds" = "#C8876EFF",
                   "Benzenoids" = "#2F533CFF",
                   "Phenylpropanoids and polyketides" = "#63764EFF",
                   "Organoheterocyclic compounds" = "#A87735FF",
                   "Lipids and lipid-like molecules" = "#F3C571FF",
                   "Organic acids and derivatives" = "#EEB043FF",
                   "Hydrocarbons" = "#E0D2B5FF")
# Add row anotation to HeatMap
hm_row_ann <- rowAnnotation(`Superclass` = hm_fdata$classyfireR_Superclass,
                            `Tukey test 2` = anno_text(hm_fdata$Late,
                                                       location = 0.5,
                                                       just = "center",
                                                       gp = gpar(fill = "#F8766D",
                                                                 col = "#666666"),
                                                       show_name = TRUE),
                            `Tukey test 1` = anno_text(hm_fdata$Medium,
                                                       location = 0.5,
                                                       just = "center",
                                                       gp = gpar(fill = "#00B0F6",
                                                                 col = "#666666"),
                                                       show_name = TRUE),
                            `Tukey test 0` = anno_text(hm_fdata$Early,
                                                       location = 0.5,
                                                       just = "center",
                                                       gp = gpar(fill = "#E76BF3",
                                                                 col = "#666666"),
                                                       show_name = TRUE),
                            col = list(`Superclass` = cols_metclass),
                            show_annotation_name = T,
                            show_legend = F)
# Factor levels color
cols_levels <- c("Shade" = "#32373AFF",
                 "Early" = "#E76BF3",
                 "Light" = "#F5BC5CFF",
                 "Late" = "#F8766D",
                 "Talag" = "#4457A5FF",
                 "Medium" = "#00B0F6",
                 "Alto Tena" = "#4F9D4EFF",
                 "Alto Pano" = "#B24422FF")
cols_light <- c("Light" = "#F5BC5CFF",
                "Shade" = "#32373AFF")
cols_age <- c("Early" = "#E76BF3",
              "Medium" = "#00B0F6",
              "Late" = "#F8766D")
cols_location <- c("Alto Pano" = "#B24422FF",
                   "Alto Tena" = "#4F9D4EFF",
                   "Talag" = "#4457A5FF")
# Add top anotation to HeatMap
top_info_ann <- HeatmapAnnotation(`Age factor` = hm_pdata$Group,
                                  col = list(`Age factor` = cols_levels),
                                  show_annotation_name = T,
                                  show_legend=F, 
                                  border = TRUE)
# Color scale
mycol <- colorRamp2(c(-2, 0, 2),
                    c("blue", "white", "red"))
# Heatmap matrix plotting
hm_plot <- Heatmap(hm_height,
                   col = mycol,
                   border_gp = grid::gpar(col = "black", lty = 0.05),
                   rect_gp = grid::gpar(col = "black", lwd = 0.75),
                   clustering_distance_columns = "euclidean",
                   clustering_method_columns = "complete",
                   top_annotation = top_info_ann,
                   right_annotation = hm_row_ann,
                   row_names_max_width = unit(10, "cm"),
                   show_heatmap_legend = FALSE,
                   row_km = 3, column_km = 2,
                   row_title = c("a", "b", "c"))
# Color scale legend
lgd1 <- Legend(col_fun = mycol,
               title = "glog abundance",
               direction = "horizontal" )
# Factor levels legend
lgd2 <- Legend(labels = c("Light (Pos)",
                          "Shade (Neg)"),
               legend_gp = gpar(fill = cols_light),
               title = "Light factor", ncol = 1)
lgd3 <- Legend(labels = c("Early (0)",
                          "Medium (1)",
                          "Late (2)"),
               legend_gp = gpar(fill = cols_age),
               title = "Age factor", ncol = 1)
lgd4 <- Legend(labels = c("Alto Pano (B)",
                          "Alto Tena (C)",
                          "Talag (A)"),
               legend_gp = gpar(fill = cols_location),
               title = "Location factor", ncol = 1)
# Metabolite class Legend
lgd5 <- Legend(labels = c(unique(hm_fdata$classyfireR_Superclass)),
               legend_gp = gpar(fill = cols_metclass),
               title = "Metabolite superclass", ncol = 3)
# Converting to ggplot
gg_heatmap <- grid.grabExpr(draw(hm_plot))
gg_heatmap <- ggpubr::as_ggplot(gg_heatmap)
# Legends
all_legends <- packLegend(lgd1,
                          #lgd2,
                          lgd3,
                          #lgd4,
                          lgd5,
                          direction = "horizontal")
gg_legend <- grid.grabExpr(draw(all_legends))
gg_legend_fn <- ggpubr::as_ggplot(gg_legend)
# Heatmap plot
gcms_hm <- plot_grid(gg_legend_fn,
                     gg_heatmap, ncol = 1,
                     rel_heights = c(0.045, 0.880))
gcms_hm
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
#ggsave(filename = "Result/notame_Result/figure_2_to_JFCA.pdf", plot = figure_two,
#       width = 175, height = 155, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
#ggsave(filename = "Result/notame_Result/figure_2_to_JFCA.png", plot = figure_two,
#       width = 175, height = 155, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.jpg) file
#ggsave(filename = "Result/notame_Result/figure_2_to_JFCA.jpg", plot = figure_two,
#       width = 175, height = 155, units = "mm", dpi = 300, scale = 2.5)




# Removal sample group from dataset
volc_pqn <- drop_flagged(pqn_set)[, drop_flagged(pqn_set)$Analysis_group != "Sample"]

# Drop QC
volc_pqn <- drop_qcs(volc_pqn)

# Extracting SubQC light factor
volc_pqn <- volc_pqn[, volc_pqn$Light_factor != "Other factor"]
pData(volc_pqn) <- droplevels(pData(volc_pqn))

# Fold change between light and shade levels
volc_fc <- fold_change(volc_pqn, group = "Light_factor")

# two-sample t-test performing
volc_t <- perform_t_test(volc_pqn, formula_char = "Feature ~ Light_factor")

# Adding the fold change to the t-test data
volc_data <- left_join(volc_t, volc_fc)

## Log-transform for visualization
volc_data$logP <- -log10(volc_data$Light_vs_Shade_t_test_P)
volc_data$log2_fc <- log2(volc_data$Shade_vs_Light_FC)

# Determine point colors based on significance and fold change
volc_data <- volc_data %>%
  mutate(point_color = case_when(
    Light_vs_Shade_t_test_P < 0.05 & log2_fc < -1 ~ "down", # significantly down
    Light_vs_Shade_t_test_P < 0.05 & log2_fc > 1 ~ "up",   # significantly up
    log2_fc < 1 & log2_fc > -1 ~ "No fc", # fold change < 2
    Light_vs_Shade_t_test_P > 0.05 & log2_fc < -1 | Light_vs_Shade_t_test_P > 0.05 & log2_fc > 1 ~ "fc" # not significant and fold change > 2
  ))

# Creating a new small table of the annotated compounds
volc_compouds <- subset(volc_data, Light_vs_Shade_t_test_P < 0.05) # Keep metabolites with p-value < 0.05
volc_compouds <- left_join(meta_table, volc_compouds)

# Volcano plot
  ggplot(volc_data, aes(log2_fc, logP, color = point_color)) +
    geom_point(size = 2,
               alpha = 0.4) +
    scale_colour_manual(values = c("No fc" = "#808080",
                                   "down" = "#8080ff",
                                   "up" = "#ff80ff",
                                   "fc" = "#ff8080"
    )) +
    theme_classic() +
    theme(legend.position = "none") +
    #geom_point(data = volc_compouds,
    #           aes(shape = meta_table$Identification_level,
    #               color = meta_table$Identification_level),
    #           size = 2) +
    #labs(shape = 'Identification level',
    #     color = 'Identification level') +
    ggrepel::geom_label_repel(data = volc_compouds,
                              aes(label = meta_table$Metabolite),
                              color = "black",
                              box.padding = 0.37,
                              label.padding = 0.22,
                              label.r = 0.30,
                              cex = 2.5,
                              max.overlaps = 20) +
    xlab(bquote(Log[2](Shade/Light))) +
    ylab(bquote(-Log[10](p-value))) +
    scale_y_continuous(limits = c(0,4), labels = function(i) 10^-i) +
    #theme(legend.position = c(0.945, 0.913),
    #      legend.background = element_rect(fill = "white", color = "black")) +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(fill= "transparent")) +
    geom_vline(xintercept = 1, linetype = "longdash", colour="gray") +
    geom_vline(xintercept = -1, linetype = "longdash", colour="gray") +
    geom_hline(yintercept = -log10(0.05), linetype = "longdash", colour="gray")

# Save plot
#ggsave('Result/notame_Result/figure_s1b.pdf',
#       width = 5, height = 4, device='pdf', dpi="print")



######################### SAMPLE ################




  # Removal SubQC group from dataset
  volc_pqn <- drop_flagged(pqn_set)[, drop_flagged(pqn_set)$Analysis_group != "Sub_Quality_Control"]
  
  # Drop QC
  volc_pqn <- drop_qcs(volc_pqn)
  
  # Extracting sample light factor
  #volc_pqn <- volc_pqn[, volc_pqn$Light_factor != "Other factor"]
  #pData(volc_pqn) <- droplevels(pData(volc_pqn))
  
  # Fold change between light and shade levels
  volc_fc <- fold_change(volc_pqn, group = "Light_factor")
  
  # two-sample t-test performing
  volc_t <- perform_t_test(volc_pqn, formula_char = "Feature ~ Light_factor")
  
  # Adding the fold change to the t-test data
  volc_data <- left_join(volc_t, volc_fc)
  
  ## Log-transform for visualization
  volc_data$logP <- -log10(volc_data$Light_vs_Shade_t_test_P)
  volc_data$log2_fc <- log2(volc_data$Shade_vs_Light_FC)
  
  # Determine point colors based on significance and fold change
  volc_data <- volc_data %>%
    mutate(point_color = case_when(
      Light_vs_Shade_t_test_P < 0.05 & log2_fc < -1 ~ "down", # significantly down
      Light_vs_Shade_t_test_P < 0.05 & log2_fc > 1 ~ "up",   # significantly up
      log2_fc < 1 & log2_fc > -1 ~ "No fc", # fold change < 2
      Light_vs_Shade_t_test_P > 0.05 & log2_fc < -1 | Light_vs_Shade_t_test_P > 0.05 & log2_fc > 1 ~ "fc" # not significant and fold change > 2
    ))
  
  # Creating a new small table of the annotated compounds
  volc_compouds <- subset(volc_data, Light_vs_Shade_t_test_P < 0.05) # Keep metabolites with p-value < 0.05
  volc_compouds <- left_join(meta_table, volc_compouds)
  
  # Volcano plot
  ggplot(volc_data, aes(log2_fc, logP, color = point_color)) +
    geom_point(size = 2,
               alpha = 0.4) +
    scale_colour_manual(values = c("No fc" = "#808080",
                                   "down" = "#8080ff",
                                   "up" = "#ff80ff",
                                   "fc" = "#ff8080"
    )) +
    theme_classic() +
    theme(legend.position = "none") +
    #geom_point(data = volc_compouds,
    #           aes(shape = meta_table$Identification_level,
    #               color = meta_table$Identification_level),
    #           size = 2) +
    #labs(shape = 'Identification level',
    #     color = 'Identification level') +
    ggrepel::geom_label_repel(data = volc_compouds,
                              aes(label = meta_table$Metabolite),
                              color = "black",
                              box.padding = 0.37,
                              label.padding = 0.22,
                              label.r = 0.30,
                              cex = 2.5,
                              max.overlaps = 20) +
    xlab(bquote(Shade/Light (Log[2]))) +
    ylab(bquote(P-value (Log[10]))) +
    scale_y_continuous(limits = c(0,4), labels = function(i) 10^-i) +
    #theme(legend.position = c(0.945, 0.913),
    #      legend.background = element_rect(fill = "white", color = "black")) +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(fill= "transparent")) +
    geom_vline(xintercept = 1, linetype = "longdash", colour="gray") +
    geom_vline(xintercept = -1, linetype = "longdash", colour="gray") +
    geom_hline(yintercept = -log10(0.05), linetype = "longdash", colour="gray")



  






  
  



# Drop QC
tk_no_qc <- drop_qcs(hm_scl_set)
# Drop samples
tk_no_sample <- tk_no_qc[, tk_no_qc$Analysis_group != "Sample"]
pData(tk_no_sample) <- droplevels(pData(tk_no_sample))
# Extracting levels of age factor
tk_age <- tk_no_sample[, tk_no_sample$Age_factor != "Other factor"]
pData(tk_age) <- droplevels(pData(tk_age))
# Extracting identified metabolites
tk_age_set <- tk_age[!is.na(tk_age@featureData@data$Metabolite),]



# Extracting feature height table
tk_height <- exprs(tk_age_set)
# Extracting sample information
tk_pdata <- tk_age_set@phenoData@data
# Extracting feature information
tk_fdata <- tk_age_set@featureData@data
# Adding row and column name
#rownames(tk_height) <- tk_fdata$Metabolite
colnames(tk_height) <- tk_pdata$Group


metabo1 <- t(tk_height)



tk_age_var <- colnames(tk_height)

metabo1 <- data.frame(tk_age_var, metabo1)





# Bartlett test of homogeneity of variances
# Batch test of all annotated features
# create two empty vectors to add values to
Feature_ID <- c()
bartlett_pValue <- c()
# write for loop
for(i in 2:79) {
  bar_res <- bartlett.test(metabo1[, i] ~ metabo1$tk_age_var)    # assign for loop iteration to variable
  bartlett_pValue <- c(bartlett_pValue, bar_res[["p.value"]])     # append current value to vector
  Feature_ID <- c(Feature_ID, names(metabo1[i]))     # append to vector
  tk_bar_res <- data.frame(Feature_ID, bartlett_pValue)     # build dataframe
}


tk_bar_res <- tk_bar_res %>% mutate(bartlett_group = ifelse(bartlett_pValue > 0.05, "equal variances", "Unequal variances"))


# ANOVA test
# Batch test of all annotated features
# create two empty vectors to add values
Feature_ID <- c()
anova_pValue <- c()
# write for loop
for(i in 2:79) {
  anova_model <- glm(metabo1[, i] ~ metabo1$tk_age_var)
  anova_res <- anova(anova_model, test="LRT")                  # assign for loop iteration to variable
  anova_pValue <- c(anova_pValue, anova_res$`Pr(>Chi)`[2])                     # append current value to vector
  Feature_ID <- c(Feature_ID, names(metabo1[i]))               # append to vector
  tk_anova_res <- data.frame(Feature_ID, anova_pValue)           # build dataframe
}


tk_anova_res <- tk_anova_res %>% mutate(anova_group = ifelse(anova_pValue > 0.05, "Not Significant", "Significant"))




# agricolae package installation and library loadding
#install.packages("agricolae", repos = "https://cran.r-project.org")
library(agricolae)




# Tukey test
# Batch test of all annotated features
# create two empty vectors to add values
Feature_ID <- c()
tk_group_early <- c()
tk_group_late <- c()
tk_group_medium <- c()
# write for loop
for(i in 2:79) {
  tk_model <- aov(metabo1[, i] ~ metabo1$tk_age_var)
  tk_res <- HSD.test(tk_model, "metabo1$tk_age_var", group=TRUE, console=TRUE)   # assign for loop iteration to variable
  tk_group_early <- c(tk_group_early, tk_res[["groups"]]["Early",2])
  tk_group_late <- c(tk_group_late, tk_res[["groups"]]["Late",2])
  tk_group_medium <- c(tk_group_medium, tk_res[["groups"]]["Medium",2])
  Feature_ID <- c(Feature_ID, names(metabo1[i]))               # append to vector
  tk_tk_res <- data.frame(Feature_ID, tk_group_early, tk_group_late, tk_group_medium)           # build dataframe
}

# Summary of Tukey test result
tk_result <- left_join(tk_bar_res, tk_anova_res)
tk_result <- left_join(tk_result, tk_tk_res)
#rownames(tk_result) <- tk_result$Feature_ID
#tk_result$row_name <- rownames(tk_result)

# Delete compounds with unequal variances
tk_result_filter <- tk_result[tk_result$bartlett_pValue > 0.05, ]
# Delete ANOVA Not Significant compounds
tk_result_filter <- tk_result_filter[tk_result_filter$anova_pValue < 0.05, ]
# Delete compounds with equal tukey group
tk_result_filter$tk_summ <- paste(tk_result_filter$tk_group_early, tk_result_filter$tk_group_late, tk_result_filter$tk_group_medium, sep = "")
tk_result_filter <- tk_result_filter[tk_result_filter$tk_summ != "aaa", ]
# Adding all annotated features
tk_result = subset(tk_result, select = -c(tk_group_early, tk_group_late, tk_group_medium, bartlett_pValue,
                                          bartlett_group, anova_pValue, anova_group) )
tk_result_filter <- right_join(tk_result_filter, tk_result)








# Table to heat map
hm_tk_table <- data.frame(Feature_ID = tk_result_filter$Feature_ID,
                          Early = tk_result_filter$tk_group_early,
                          Late = tk_result_filter$tk_group_late,
                          Medium = tk_result_filter$tk_group_medium
                          #Early = tk_result_filter$tk_group_early,
                          #Late = tk_result_filter$tk_group_late,
                          #Medium = tk_result_filter$tk_group_medium,
                          #Early = tk_result_filter$tk_group_early,
                          #Late = tk_result_filter$tk_group_late,
                          #Medium = tk_result_filter$tk_group_medium
                          )
hm_tk_table[is.na(hm_tk_table)]  <- "ns"


#rownames(hm_tk_table) <- hm_tk_table$Row_name


#hm_tk_table <- as.matrix(hm_tk_table)
#hm_tk_table <- hm_tk_table[,-1]


#colnames(hm_tk_table)  <- c("Early", "Late", "Medium", "Early", "Late",
#                            "Medium", "Early", "Late", "Medium")
#hm_tk_table[is.na(hm_tk_table)]  <- "ns"








metabo1 <- t(tk_height)


tk_age_var <- colnames(tk_height)

# Data extraction
metabo1 <- tk_height["Coniferyl alcohol",]

metabo1 <- data.frame(tk_age_var, metabo1)

#install.packages("tidyverse", type="source")
library(tidyverse)

metabo1 <- metabo1 %>%
  pivot_longer(-tk_age_var, names_to = "variables", values_to = "value")


metabo1 <- metabo1 %>%
  group_by(variables) %>%
  bartlett.test(value ~ tk_age_var, metabo1)



stat.test <- mydata.long %>%
  group_by(variables) %>%
  t_test(value ~ Species) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test



# Boxplot
ggplot(metabo1, aes(x=colnames(tk_height),
                    y=metabo1, fill=colnames(tk_height))) +
  geom_boxplot() +
  ggtitle("Coniferyl alcohol") +
  guides(x=guide_axis(title = NULL),
         y=guide_axis(title = "glog(Peak height)"),
         fill=guide_legend(title="Age factor"))



bar_res <- bartlett.test(metabo1 ~ tk_age_var, data = metabo1)





    
  
  
  

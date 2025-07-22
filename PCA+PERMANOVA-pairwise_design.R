rm(list=ls()) #clear workspace

library("readxl")
library("dplyr")
library("VIM") # KNN Imputation
library ("ggplot2")
library ("ggcorrplot") # correlation plot 
library ("factoextra") # PCA visualization
library ("vegan") # PERMANOVA

# LOAD DATA #
# here: mean per tree dataset
field_import <- read_excel("./YourDATA.xlsx", na = "na")

# PREP DATA #
# KNN Imputation, to get values for NAs 
field_imputed <- kNN(field_import, k = 10)
# To remove the '.imp' columns created by kNN
field_imputed <- field_imputed %>%
  select(-ends_with("_imp"))

# only one site 
field_imputed <- field_imputed %>%
  filter(site == "Site1")

# rename columns 
field_imputed <- field_imputed %>%
  rename("AC" = "air_capacity",
         "BD" = "bulk_density", 
         "MacroP - 10" = "macroporosity_cores", 
         "Pore dia - 10" = "pore_diameter_cores", 
         "POM - 10" = "POM_cores",
         "MacroP - 3" = "macroporosity_ss",
         "Pore dia - 3" = "pore_diameter_ss", 
         "POM - 3 " = "POM_ss",
         "pH" = "pH",
         "CN" = "CN",
         "PA" = "total_PA", 
         "Dehydrogenase" = "dehydrogenase", 
         "Shannon" = "shannon"
  )

# merge site+growth+treatment in one column 
field_imputed$interaction <- interaction (field_imputed$site, field_imputed$growth, field_imputed$treatment)

# omit non-numeric columns, redundant columns 
field_edited <- field_imputed [,-c(1:5, 7, 18, 19, 20, 21, 25)]

# sava data as matrix 
field <- as.matrix(field_edited)

# normalize the data 
field_nor <- scale (field)

# CORRELATION MATRIX #
corr_matrix <- cor (field_nor)

# correlation plot 
corr_plot <- ggcorrplot(
  corr_matrix,
  type = "lower",       
  lab = TRUE,           
  lab_size = 3,          
  outline.color = "gray"
) +
  scale_fill_gradientn(
    colors = c("darkblue", "lightblue", "white", "lightcoral", "darkred"),
    values = scales::rescale(c(-1, -0.6, 0, 0.6, 1)),
    limits = c(-1, 1)
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),             
    axis.text.x = element_text(angle = 45,
                               hjust = 1, 
                               vjust = 1)
  )
corr_plot

ggsave(filename = "Site1_corr_matrix.png", plot = corr_plot, width = 8, height = 6, dpi = 300)

# COMPUTE PCA #
field.pca <- princomp (field_nor)
summary (field.pca)

# visualize as screeplot
fviz_eig (field.pca, addlabels = TRUE)

# PERMANOVA #
# to add the metadata columns (site, treatment, block..) to the data 
selected_columns <- field_imputed %>% select(1:5)
field_meta_edited <- field_edited %>% bind_cols(selected_columns)

# make factors out of character vectors 
field_meta_edited <- field_meta_edited %>% 
  as_tibble() %>%
  mutate(across(where(is.character), as.factor))

# compute distance matrix 
perm_dist <- vegdist(field_edited, method='euclidean')

# check for homogeneity of multivariate dispersion as a prerequisite for PERMANOVA (p>0.05 necessary for application)
# type = median is set here because it is recommended by the package authors
dispersion <- betadisper(perm_dist, group = field_imputed$interaction, type = "median")
anova (dispersion)

# 1. whole plot PERMANOVA - TREATMENT EFFECT # 
# 1.1 define permutation object 
permut_whole <- how (within = Within (type = "none"), 
                     plots = Plots (strata = field_meta_edited$pair,type = "free"),
                     nperm = 3999, 
                     observed = TRUE)
# keeping pairs intact, to simulate treatment effects between pairs
# permutations is here set to 3999, as testing showed this provides relatively stable results

# 1.2 create a set of permutations 
perms <- rbind(1:nrow(field_meta_edited),
               shuffleSet(n = nrow(field_meta_edited), control = permut_whole, nset = 3999))
perms
# to get many possible arrangements and with that a robust stats outcome

# 1.3 create object to store results in 
results <- matrix(nrow = nrow(perms), ncol = 4)
colnames(results) <- c("treatment", "pair", "residual", "total")

# 1.4 loop through permutations  
for (i in 1:nrow(perms)) {
  temp.data <- field_meta_edited[perms[i, ], ]
  temp <- adonis2(perm_dist ~ treatment + pair,
                  data = temp.data,
                  by = "terms",
                  method = "euclidean",
                  permutations = 0)
  results[i, ] <- t(temp$SumOfSqs)
}
head (results)
# permutations are here set to zero because they are manually handled 

# 1.5 pseudo-F-statistics
results <- results |>
  data.frame() |>
  mutate(df.treatment = temp$Df[1],
         df.pair = temp$Df[2]) |>
  mutate(
    F.treatment = (treatment / df.treatment) / (pair / df.pair),  
    R2 = treatment / total 
  )
results[1,]

# 1.6 acquire p value
pvalue_tr <- with(results, round(sum(F.treatment >= F.treatment[1]) / length(F.treatment), digits = 3))
pvalue_tr

# 2. split-plot PERMANOVA - GROWTH EFFECT # 
# 2.1 define permutation object
permut_split <- how(within = Within(type = "free"),
                    plots = Plots(type = "none"),
                    blocks = as.factor(field_meta_edited$pair),
                    observed = TRUE,
                    nperm = 1024)
# within-pair permutations while keeping pairs separate 
# 1024 is here the maximum number of permutations (10^2) 

# 2.2 compute split-plot PERMANOVA
perm_result_res <- adonis2(perm_dist ~ pair + growth,
                           data = field_meta_edited,
                           permutations = permut_split
)
perm_result_res

# PLOT #
# save PERMANOVA output as objects 
r_squared_gr <- perm_result_res$R2[2]
pvalue_gr <- perm_result_res$`Pr(>F)`[2]
r_squared_tr <- results$R2[1]

# to plot the pair-lines
pca_coords <- get_pca_ind(field.pca)$coord  
pca_coords <- as.data.frame(pca_coords)
pca_coords$sample <- rownames(pca_coords)  

# matching identifier column for "field_imputed"
field_imputed$sample <- 1:20
field_export <- merge(field_imputed, pca_coords, by = "sample")

# convert to factor, to plot lines 
field_export$pair<- as.factor(field_export$pair)

# Score plot #
scores <- fviz_pca_biplot (field.pca, habillage = field_imputed$interaction, label = "none", addEllipses = TRUE,  ellipse.type = "confidence", point = 16, pointsize = 4, mean.point = FALSE, labelsize = 7, repel = TRUE, col.var = NA, alpha.ind = 1, ellipse.alpha = 0.3, arrowsize = 0.75) +
  theme_light (base_size = 25) +
  theme (panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         plot.title = element_blank(), 
         legend.title = element_blank(),
         legend.position ="none") +
  scale_color_manual(values = c("goldenrod3", "orangered4", "steelblue4", "honeydew4") 
  ) +
  scale_fill_manual(values = c("goldenrod3", "orangered4", "steelblue4", "honeydew4"),
                    labels = c("Control - Better", "Control - Worse", "Tagetes - Better", "Tagetes - Worse")) +
  guides (color = "none",
          fill = guide_legend(override.aes = list (size = 3, alpha = 1))) +
  annotate("text", x = -3.5, y = 5, label = expression(italic("Treatment")), size = 8.5,
           color = "black", hjust = 0) + 
  annotate("text", x = 1.1, y = 5, label = expression(italic("Growth")), size = 8.5,
           color = "black", hjust = 0) +
  annotate("text", x = -4.2, y = 4.3, label = expression("R² ="), size = 9,
           color = "black", hjust = 0) + 
  annotate("text", x = -4.2, y = 3.6, label = expression("p ="), size = 9,
           color = "black", hjust = 0) +
  annotate("text", x = 0.3, y = 4.3, label = expression("R² ="), size = 9,
           color = "black", hjust = 0) + 
  annotate("text", x = 0.3, y = 3.6, label = expression("p ="), size = 9,
           color = "black", hjust = 0) +
  annotate("text", x = -3.2, y = 4.3, label = paste(round(r_squared_tr, 3)), size = 9,
           color = "black", hjust = 0) +
  annotate("text", x = -3.2, y = 3.6, label = paste(round(pvalue_tr, 3)), size = 9,
           color = "black", hjust = 0) +
  annotate("text", x = 1.2, y = 4.3, label = paste(round(r_squared_gr, 3)), size = 9,
           color = "black", hjust = 0) +
  annotate("text", x = 1.2, y = 3.6, label = paste(round(pvalue_gr, 3)), size = 9,
           color = "black", hjust = 0) +
  geom_line (data = field_export, aes (x = Dim.1, y = Dim.2, group = pair)) +
  coord_fixed (ratio = 0.7)
scores

# to increase line thickness of ellipses
q <- ggplot_build(scores)
q$data[[2]]$linewidth = 1.3
q <- ggplot_gtable(q)
plot(q)

ggsave("Site1_scores.png", width= 30, height = 17.5, units = "cm")

# Loadings Plot #
# create a grouping factor for variables
var_groups <- factor(c("retention", "physico-chemical", "structure", 
                       "structure", "structure", "structure", "structure", "structure", "physico-chemical", "physico-chemical", "physico-chemical", "PA", "MO", "MO"))
# custom color palette
custom_colors <- c("retention" = "darkgrey",  
                   "structure" = "hotpink3",
                   "physico-chemical" = "blue3", 
                   "PA" = "darkorange1", 
                   "MO" = "olivedrab4")

# visualize variables
loadings <- fviz_pca_var(field.pca, col.var = var_groups, repel = TRUE, labelsize = 8) +
  scale_color_manual(values = custom_colors) +
  theme_light (base_size = 25) +
  theme(
    plot.title = element_blank(),      
    panel.grid = element_blank(),      
    axis.text = element_text(size = 18), 
    axis.title = element_text(size = 24), 
    legend.title = element_blank(), 
    legend.position ="none")
loadings

ggsave(filename = "Site1_loadings.png", plot = loadings, bg = "white", width = 8.85, height = 6.85, dpi = 300)



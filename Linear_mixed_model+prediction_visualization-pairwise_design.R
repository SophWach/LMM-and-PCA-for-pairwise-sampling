rm(list=ls()) # clear workspace

library("readxl")
library("ggplot2")
library ("dplyr")
library ("lme4") # to fit linear mixed model
library("emmeans") # to acquire p-values
library ("performance") # to check model quality and assumptions
library("DHARMa") # for simulation-based diagnostics 
library ("ggeffects") # to compute model predictions 
library ("ggpubr") # to add brackets to ggplot
library ("cowplot") # to draw on the ggplot

# LOAD DATA #
# here: mean per tree dataset
model_data <- read_excel("./YourData.xlsx")

# filter for one site and define growth and treatment as factors 
model_data_Site1 <- 
  model_data %>%
  filter (site == "Site1") %>%
  mutate (growth = as.factor (growth),
          treatment = as.factor (treatment), 
          pair = as.factor (pair))
model_data_Site1

# get tree group columns (merging growth and treatment columns)
model_data_Site1$treatmentgrowth <- interaction (model_data_Site1$growth, model_data_Site1$treatment)

# FIT LINEAR MIXED MODEL #
# treatment and growth as interacting fixed effects, allowing for a unique intercept per pair 
model.Variable1 <- lmer (Variable1 ~ treatment * growth + (1 | pair), data = model_data_Site1)
summary (model.Variable1)

# ACQUIRE P VALUES
# in-treatment
# compute estimated marginal means for each combination of treatment and growth for Variable1
emm <- emmeans(model.Variable1, ~ treatment * growth)
# perform pairwise comparisons among all treatment-growth combinations
pairs(emm)

# between treatments 
# compute estimated marginal means for the treatment factor of Variable1, averaging over growth
emm_treatment <- emmeans(model.Variable1, ~ treatment)
# perform pairwise comparison between the two treatments 
contrast(emm_treatment, method = "pairwise")

# DIAGNOSTICS # 
# check model for normality,linearity, homogeneity 
check_model(model.Variable1, check = c("qq", "linearity", "homogeneity"), detrend = FALSE)

ggsave("Site1_Diagnostics_Variable1.png", width= 30, height = 17.5, units = "cm")

# additionally simulation-based residual diagnostics
simulateResiduals(model.Variable1) %>% plot
simulateResiduals(model.Variable1) %>% plot(rank = FALSE)

# EFFECT PLOTS #
# to extract the modes of the random effects; so how each pairs intercepts varies from the models overall intercept 
random_effects <- ranef(model.Variable1)
random_effects

# 1. RANDOM INTERCEPTS #
random_effects_df <- as.data.frame(random_effects$pair)
random_effects_df$pair <- factor(rownames(random_effects_df), levels = 1:10) # levels is here set to 10, because of 10 tree pairs/site
confint <- sqrt(attr(random_effects$pair, "postVar")[1, , ])  # Random effects confidence intervals

ggplot(random_effects_df, aes(x = pair, y = `(Intercept)`)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = `(Intercept)` - 1.96 * confint, 
                    ymax = `(Intercept)` + 1.96 * confint), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Random Intercepts with Confidence Intervals",
       x = "Pair",
       y = "Random Intercept (Deviation from Overall Intercept)") +
  theme_minimal()

ggsave("Site1_Random_intercepts_Variable1.png", width= 20, height = 10, units = "cm")

# 2. PREDICTIONS FOR THE RESPONSE VARIABLE # 
# compute model predictions per tree group (estimated marginal means + confidence intervals)
pr <- predict_response(model.Variable1, terms = c("treatment", "growth"), type = "random", interval = "confidence")
as.data.frame(pr)
# get tree group column 
pr$treatmentgrowth <- interaction (pr$group, pr$x)
# rename column that contains predicted means (for plotting)
pr <- pr %>%
  rename(Variable1 = predicted)

# Mean-Jitter Plot #
plot_stats <- ggplot(model_data_Site1, aes(treatmentgrowth, Variable1, color = treatmentgrowth)) +
  geom_point(data = pr, size = 6.5) + # mean points = estimated marginal mean 
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.9, dodge.width = 0.8), alpha = 0.5, size = 5) + # jitters = measured values per tree
  geom_errorbar(
    aes(ymin = pmax(0,conf.low), ymax = conf.high), data = pr, 
    width = 0.2, position = position_dodge(0.8), linewidth = 1.3
  ) + # error bars = 95% confidence interval of estimated marginal mean
  scale_color_manual(values = c("goldenrod3", "orangered4", "steelblue4", "honeydew4"), 
                     labels = c("Control - Better", "Control - Worse", "Tagetes - Better", "Tagetes - Worse")) +
  theme_light(base_size = 28) +
  labs(y = "Variable1") + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0)) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 28),
        axis.text.y = element_text(size = 28),
        legend.position = "none") +
  scale_x_discrete(labels = c("better.control" = "\nBetter",
                              "worse.control" = "\nWorse",
                              "better.tagetes" = "\nBetter",
                              "worse.tagetes" = "\nWorse"))

# add statistic results to the plot #

# create an object that can be rendered
q <- ggplot_build(plot_stats)
q$data[[1]]$x[1:2]

# Calculate midpoints between our 4 elements
midpoint1 <- mean(q$data[[1]]$x[1:2])
midpoint2 <- mean(q$data[[1]]$x[3:4])

# calculate extended xmin and xmax for the larger bracket
extension_left <- 0.2 * (q$data[[1]]$x[2] - q$data[[1]]$x[1]) 
extension_right <- 0.2 * (q$data[[1]]$x[4] - q$data[[1]]$x[3]) 

# add brackets to the plot
plot_stats_with_brackets <- plot_stats +
  # first bracket
  geom_bracket(
    xmin = q$data[[1]]$x[1], xmax = q$data[[1]]$x[2],  
    y.position = 2, label.size = 9, label = "0.002", # label is here the controlbetter-controlworse p value generated by emmeans () in line 39
    inherit.aes = FALSE, tip.length = 0.01, size = 0.7
  ) +
  # second bracket
  geom_bracket(
    xmin = q$data[[1]]$x[3], xmax = q$data[[1]]$x[4],  
    y.position = 2, label.size = 9, label = "0.87",  # label is here the tagetesbetter-tagetesworse p value generated by emmeans () in line 39
    inherit.aes = FALSE, tip.length = 0.01, size = 0.7
  ) +
  # third bracket
  geom_bracket(
    xmin = midpoint1 - extension_left,  
    xmax = midpoint2 + extension_right,  
    y.position = 2.25, label.size = 9, label = "0.002", # label is here the p value generated by emmeans () in line 45
    inherit.aes = FALSE, tip.length = 0.01, size = 0.7
  )

# add tagetes and control annotation 
final_plot <- ggdraw(plot_stats_with_brackets) +
  draw_label("Control", x = 0.375, y = 0.12, size = 28) +
  draw_label("Tagetes", x = 0.764, y = 0.12, size = 28)
final_plot

ggsave("Site1_mean_plot_Variable1.png", width= 20, height = 17.5, units = "cm")


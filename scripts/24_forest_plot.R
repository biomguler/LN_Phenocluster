#############################################################################
#Title:    Forest plots
#Function: Create custom forest plots
#Author:   Murat Guler
#Time:     May 17th 2024
#############################################################################
#source: https://www.khstats.com/blog/forest-plots/#just-the-code
#############################################################################
rm(list=ls())
gc()
#############################################################################
#required packages
library(data.table)
library(tidyverse)
library(patchwork)
library(gt)

data <- read.table(text = "
rsID NEA EA OR CI_95_low CI_95_up P_value Study
rs13018297 G C 2.19 1.68 2.86 9.36E-09 'Asset-1sided'
rs13018297 G C 1.67 1.38 2.02 1.12E-07 'LN'
rs13018297 G C 1.79 1.41 2.26 1.23E-06 'Cell-B'
rs13018297 G C 1.65 1.35 2.03 1.78E-06 'Drug-G1'
rs13018297 G C 1.94 1.47 2.56 2.53E-06 'Soma-G2  Meta analysis'
rs13018297 G C 2.63 1.65 4.18 4.64E-05 'Soma-G1'
rs13018297 G C 2.93 1.73 4.97 6.11E-05 'Soma-G1 Meta analysis'
rs13018297 G C 1.62 1.28 2.05 6.29E-05 'Soma-G2'
rs13018297 G C 2.01 1.38 2.94 2.90E-04 'DLBCL'
rs13018297 G C 1.83 1.30 2.60 6.32E-04 'Cell-P Meta analysis'
rs13018297 G C 1.80 1.28 2.53 7.78E-04 'Cell-P'
rs13018297 G C 2.69 1.45 5.00 1.79E-03 'FL'
rs13018297 G C 1.86 1.24 2.78 2.58E-03 'MM'
rs13018297 G C 3.74 1.47 9.54 5.79E-03 'LPL-WM'
rs13018297 G C 3.64 1.37 9.67 9.61E-03 'MZL'
rs13018297 G C 1.64 1.13 2.38 9.97E-03 'MM-MGUS  Meta analysis'
rs13018297 G C 1.58 1.10 2.27 1.41E-02 'MM-MGUS'
rs13018297 G C 0.73 0.26 2.01 5.39E-01 'MGUS'
rs13018297 G C 0.91 0.59 1.39 6.55E-01 'CLL'
rs13018297 G C 0.81 0.28 2.31 6.93E-01 'HL'
rs561366 G A 0.86 0.81 0.90 1.58E-09 'LN'
rs561366 G A 0.85 0.80 0.89 2.25E-09 'Drug-G1'
rs561366 G A 0.83 0.78 0.89 8.61E-09 'Soma-G2'
rs561366 G A 0.74 0.66 0.83 1.53E-07 'CLL'
rs561366 G A 0.85 0.80 0.90 3.53E-07 'Asset-1sided'
rs561366 G A 0.85 0.80 0.91 4.74E-07 'Cell-B'
rs561366 G A 0.85 0.78 0.93 3.49E-04 'Cell-P'
rs561366 G A 0.86 0.78 0.94 1.22E-03 'Cell-P Meta analysis'
rs561366 G A 0.90 0.83 0.97 4.09E-03 'Soma-G2  Meta analysis'
rs561366 G A 0.88 0.80 0.97 9.08E-03 'MM-MGUS'
rs561366 G A 0.88 0.80 0.97 1.07E-02 'MM-MGUS  Meta analysis'
rs561366 G A 0.74 0.58 0.95 1.92E-02 'LPL-WM'
rs561366 G A 0.89 0.80 0.99 2.90E-02 'MM'
rs561366 G A 0.91 0.82 1.00 5.98E-02 'DLBCL'
rs561366 G A 0.82 0.62 1.08 1.51E-01 'MGUS'
rs561366 G A 0.90 0.77 1.07 2.35E-01 'FL'
rs561366 G A 0.94 0.83 1.06 3.11E-01 'Soma-G1'
rs561366 G A 1.10 0.84 1.44 5.00E-01 'MZL'
rs561366 G A 0.95 0.83 1.10 5.15E-01 'Soma-G1 Meta analysis'
rs561366 G A 1.02 0.77 1.35 8.96E-01 'HL'
", header = TRUE, stringsAsFactors = FALSE)



# rs13018297

res_rs13018297 <- data %>%
  mutate(
    log.estimate = log(OR),
    log.conf.low = log(CI_95_low),
    log.conf.high = log(CI_95_up),
    estimate = OR,
    conf.low = CI_95_low,
    conf.high = CI_95_up,
    p.value = P_value
  ) %>% filter(rsID=="rs13018297") %>%
  select(-NEA, -EA, -OR, -CI_95_low, -CI_95_up, -P_value, -rsID) %>% 
  rename(model= Study) %>% 
  arrange(p.value) 

res_rs13018297$model <- factor(res_rs13018297$model, levels = res_rs13018297$model[order(res_rs13018297$p.value)])


# Create three plot and merge them
p_rs13018297 <- 
  res_rs13018297 |>
  ggplot(aes(y = fct_rev(model))) + 
  theme_classic()
p_rs13018297

p_rs13018297 <- p_rs13018297 +
  geom_point(aes(x=log.estimate), shape=15, size=3) +
  geom_linerange(aes(xmin=log.conf.low, xmax=log.conf.high)) 
p_rs13018297

p_rs13018297 <- p_rs13018297 +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Log Odd Ratio", y="")
p_rs13018297

p_rs13018297 <- p_rs13018297 +
  coord_cartesian(ylim=c(1,21), xlim=c(-2, 3))
p_rs13018297

p_rs13018297<- p_rs13018297 +
  annotate("text", x = 1, y = 21, label = "rs13018297-C risk") + 
  annotate("text", x = -1, y = 21, label = "rs13018297-C protective")
p_rs13018297

p_rs13018297_mid <- p_rs13018297 + 
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank())
p_rs13018297_mid


#second part
# wrangle results into pre-plotting table form
res_rs13018297_plot <- res_rs13018297 |>
  # round estimates and 95% CIs to 2 decimal places for journal specifications
  mutate(across(
    c(estimate, conf.low, conf.high),
    ~ str_pad(
      formatC(.x, format = "f", digits = 2),
      width = 4,  # Adjust width as needed
      pad = "0",
      side = "right"
    )
  ),
  # add an "-" between OR estimate confidence intervals
  estimate_lab = paste0(estimate, " (", conf.low, "-", conf.high, ")")) |>
  # format all p-values in scientific notation
  mutate(p.value = format(p.value, scientific = TRUE)) |>
  # add a row of data that are actually column names which will be shown on the plot in the next step
  bind_rows(
    data.frame(
      model = "Phenocluster/Subtype",
      estimate_lab = "Odd Ratio (95% CI)",
      conf.low = "",
      conf.high = "",
      p.value = "P"
    )
  ) |>
  mutate(model = fct_rev(fct_relevel(model, "Phenocluster/Subtype")))

glimpse(res_rs13018297_plot)

ordered_levels <- c("Phenocluster/Subtype", levels(res_rs13018297$model))
res_rs13018297_plot$model <- factor(res_rs13018297_plot$model, levels = ordered_levels)


p_rs13018297_left <-
  res_rs13018297_plot  |>
  ggplot(aes(y = fct_rev(model)))
p_rs13018297_left

p_rs13018297_left <- 
  p_rs13018297_left +
  geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold")
p_rs13018297_left

p_rs13018297_left <- 
  p_rs13018297_left +
  geom_text(
    aes(x = 1, label = estimate_lab),
    hjust = 0,
    fontface = ifelse(res_rs13018297_plot$estimate_lab == "Odd Ratio (95% CI)", "bold", "plain")
  )

p_rs13018297_left

p_rs13018297_left <-
  p_rs13018297_left +
  theme_void() +
  coord_cartesian(xlim = c(0, 2)) # need to be adjusted preciously

p_rs13018297_left

# right side of plot - pvalues

p_rs13018297_right <-
  res_rs13018297_plot  |>
  ggplot() +
  geom_text(
    aes(x = 0, y = fct_rev(model), label = p.value),
    hjust = 0,
    fontface = ifelse(res_rs13018297_plot$p.value == "P", "bold", "plain")
  ) +
  theme_void() 

p_rs13018297_right

layout <- c(
  area(t = 0, l = 0, b = 30, r = 6), # left plot, starts at the top of the page (0) and goes 30 units down and 7 units to the right
  area(t = 1, l = 7, b = 30, r = 12), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
  area(t = 0, l = 13, b = 30, r = 15) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)

# final plot arrangement
tiff("rs13018297_forest.tiff", width = 16, height = 9, units = "in", res = 300)

p_rs13018297_left + p_rs13018297_mid + p_rs13018297_right + plot_layout(design = layout)

dev.off()



###############################################################################



# rs561366

res_rs561366 <- data %>%
  mutate(
    log.estimate = log(OR),
    log.conf.low = log(CI_95_low),
    log.conf.high = log(CI_95_up),
    estimate = OR,
    conf.low = CI_95_low,
    conf.high = CI_95_up,
    p.value = P_value
  ) %>% filter(rsID=="rs561366") %>%
  select(-NEA, -EA, -OR, -CI_95_low, -CI_95_up, -P_value, -rsID) %>% 
  rename(model= Study) %>% 
  arrange(p.value) 

res_rs561366$model <- factor(res_rs561366$model, levels = res_rs561366$model[order(res_rs561366$p.value)])


# Create three plot and merge them
p_rs561366 <- 
  res_rs561366 |>
  ggplot(aes(y = fct_rev(model))) + 
  theme_classic()
p_rs561366

p_rs561366 <- p_rs561366 +
  geom_point(aes(x=log.estimate), shape=15, size=3) +
  geom_linerange(aes(xmin=log.conf.low, xmax=log.conf.high)) 
p_rs561366

p_rs561366 <- p_rs561366 +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Log Odd Ratio", y="")
p_rs561366

p_rs561366 <- p_rs561366 +
  coord_cartesian(ylim=c(1,21), xlim=c(-2, 3))
p_rs561366

p_rs561366<- p_rs561366 +
  annotate("text", x = 1, y = 21, label = "rs561366-A risk") + 
  annotate("text", x = -1, y = 21, label = "rs561366-A protective")
p_rs561366

p_rs561366_mid <- p_rs561366 + 
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank())
p_rs561366_mid


#second part
# wrangle results into pre-plotting table form
res_rs561366_plot <- res_rs561366 |>
  # round estimates and 95% CIs to 2 decimal places for journal specifications
  mutate(across(
    c(estimate, conf.low, conf.high),
    ~ str_pad(
      formatC(.x, format = "f", digits = 2),
      width = 4,  # Adjust width as needed
      pad = "0",
      side = "right"
    )
  ),
  # add an "-" between OR estimate confidence intervals
  estimate_lab = paste0(estimate, " (", conf.low, "-", conf.high, ")")) |>
  # format all p-values in scientific notation
  mutate(p.value = format(p.value, scientific = TRUE)) |>
  # add a row of data that are actually column names which will be shown on the plot in the next step
  bind_rows(
    data.frame(
      model = "Phenocluster/Subtype",
      estimate_lab = "Odd Ratio (95% CI)",
      conf.low = "",
      conf.high = "",
      p.value = "P"
    )
  ) |>
  mutate(model = fct_rev(fct_relevel(model, "Phenocluster/Subtype")))

glimpse(res_rs561366_plot)

ordered_levels <- c("Phenocluster/Subtype", levels(res_rs561366$model))
res_rs561366_plot$model <- factor(res_rs561366_plot$model, levels = ordered_levels)


p_rs561366_left <-
  res_rs561366_plot  |>
  ggplot(aes(y = fct_rev(model)))
p_rs561366_left

p_rs561366_left <- 
  p_rs561366_left +
  geom_text(aes(x = 0, label = model), hjust = 0, fontface = "bold")
p_rs561366_left

p_rs561366_left <- 
  p_rs561366_left +
  geom_text(
    aes(x = 1, label = estimate_lab),
    hjust = 0,
    fontface = ifelse(res_rs561366_plot$estimate_lab == "Odd Ratio (95% CI)", "bold", "plain")
  )

p_rs561366_left

p_rs561366_left <-
  p_rs561366_left +
  theme_void() +
  coord_cartesian(xlim = c(0, 2)) # need to be adjusted preciously

p_rs561366_left

# right side of plot - pvalues

p_rs561366_right <-
  res_rs561366_plot  |>
  ggplot() +
  geom_text(
    aes(x = 0, y = fct_rev(model), label = p.value),
    hjust = 0,
    fontface = ifelse(res_rs561366_plot$p.value == "P", "bold", "plain")
  ) +
  theme_void() 

p_rs561366_right

layout <- c(
  area(t = 0, l = 0, b = 30, r = 6), # left plot, starts at the top of the page (0) and goes 30 units down and 7 units to the right
  area(t = 1, l = 7, b = 30, r = 12), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
  area(t = 0, l = 13, b = 30, r = 15) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)

# final plot arrangement
tiff("rs561366_forest.tiff", width = 16, height = 9, units = "in", res = 300)
p_rs561366_left + p_rs561366_mid + p_rs561366_right + plot_layout(design = layout)

dev.off()
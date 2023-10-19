# ============= #
#  7. Figure 2  #
# ============= #

## Set up ----
filepath <- ""
setwd(filepath)

## Load libraries ----
library(cowplot)
library(ggstance)
library(grid)
library(tidyverse)

## Figure 2a ----
### Format data ----
# List of pathogens
pathogens <- c(
  "HSV1",
  "HSV2",
  "VZV",
  "EBV",
  "CMV",
  "HHV6A",
  "HHV6B",
  "HHV7",
  "KSHV",
  "BK",
  "JC",
  "MCV",
  "HPV16",
  "HPV18",
  "T.gondii",
  "C.trachomatis",
  "H.pylori"
)

# Results data
fig2a_res <-
  # read in data
  read.csv("results/serostatus_main_MA.csv") %>%
  # filter for serostatus results only
  filter(!(pathogen %in% c("PBI", "neuro_PBI"))) %>%
  # convert to sympercents for white matter lesions
  mutate(
    beta = ifelse(outcome == "WMHV", beta * 100 , beta),
    ci_lower = ifelse(outcome == "WMHV", ci_lower * 100  , ci_lower),
    ci_upper = ifelse(outcome == "WMHV", ci_upper * 100  , ci_upper),
    # capitalise model for presentation purposes
    model = factor(str_to_sentence(model),
                   # set order for plotting
                   levels = c("Model 3", "Model 2", "Model 1")),
    # indicate significance threshold
    threshold = factor(
      case_when(
        padj < 0.05 ~ "padj<0.05",
        padj > 0.05 & p < 0.05 ~ "p<0.05",
        p > 0.05 ~ "p>0.05"
      ),
      levels = c("p>0.05",
                 "p<0.05",
                 "padj<0.05")
    ),
    # convert bacteria/protozoa to longer name
    pathogen = case_when(
      as.character(pathogen) == "Ct" ~ "C.trachomatis",
      as.character(pathogen) == "Hp" ~ "H.pylori",
      as.character(pathogen) == "Tg" ~ "T.gondii",
      TRUE ~ as.character(pathogen)
    ),
    # amend outcome name for presentation purposes
    outcome = case_when(
      outcome == "brain volume" ~ "A) Brain volume",
      outcome == "hippocampal volume" ~ "B) Hippocampal volume",
      outcome == "WMHV" ~ "C) White matter lesions"
    )
  ) %>%
  mutate(
    # set order
    pathogen = factor(pathogen, levels = pathogens),
    # indicate pathogen family
    family = case_when(
      pathogen %in% c(
        "HSV1",
        "HSV2",
        "VZV",
        "EBV",
        "CMV",
        "HHV6A",
        "HHV6B",
        "HHV7",
        "KSHV"
      ) ~ "Herpesviruses",
      pathogen %in% c("BK", "JC", "MCV") ~ "Human Polyomaviruses",
      pathogen %in% c("HPV16", "HPV18") ~ "Human Papillomaviruses",
      pathogen %in% c("T.gondii", "C.trachomatis", "H.pylori") ~ "Bacteria/Protozoa"
    )
  ) %>%
  # select cols of interest
  select(pathogen,
         family,
         outcome,
         model,
         beta,
         p,
         ci_lower,
         ci_upper,
         threshold)

### Plot - brain volume ----
fig2a_bv <- fig2a_res %>%
  # filter for brain volume outcome
  filter(outcome == "A) Brain volume") %>%
  # plot pathogen on y axis
  ggplot(aes(y = pathogen)) +
  # plot vals on x axis,
  geom_pointrangeh(
    aes(
      x = beta,
      xmin = ci_lower,
      xmax = ci_upper,
      # color by model
      col = model,
      # shape and size by threshold
      shape = threshold,
      size = threshold
    ),
    # plot models below each other
    position = position_dodgev(height = 0.75)
  ) +
  # specify number of breaks on x-axis
  scale_x_continuous(n.breaks = 6) +
  # reverse plotting of pathogens
  scale_y_discrete(limits = rev) +
  # set shape for plotting estimates (indicates threshold)
  scale_shape_manual(values = c(4, 1, 19)) +
  # set size for plotting estimate (indicates threshold)
  scale_size_manual(values = c(0.13, 0.29, 0.29)) +
  # set colors for plotting models
  scale_colour_manual(values = c("#508CC8",
                                 "#FF4A4F",
                                 "#2D506E")) +
  # include vertical dashed line at 0
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    size = 0.4,
    color = "dark grey"
  ) +
  # split into panels indicating pathogen family
  facet_grid(family ~ outcome, scales = "free", space = "free_y") +
  # amend axis labels
  labs(x = "Beta, 95% CI (ml)", y = "") +
  # set theme
  theme_classic() +
  theme(
    # set axis text size
    axis.text = element_text(size = 10),
    # amend space between facet panels
    panel.spacing = unit(0.8, "lines"),
    # format x strip text
    strip.text.x =  element_text(size = 10.5,  face = "bold" , hjust = 0.1),
    # remove y strip text & background
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    # remove legend
    legend.position = "none"
  )

### Plot - hippocampal volume ----
fig2a_hv <- fig2a_res %>%
  # filter for hippocampal volume outcome
  filter(outcome == "B) Hippocampal volume") %>%
  ggplot(aes(y = pathogen)) +
  geom_pointrangeh(
    aes(
      x = beta,
      xmin = ci_lower,
      xmax = ci_upper,
      col = model,
      shape = threshold,
      size = threshold
    ),
    position = position_dodgev(height = 0.75)
  ) +
  scale_x_continuous(n.breaks = 6) +
  scale_y_discrete(limits = rev) +
  scale_shape_manual(values = c(4, 1, 19)) +
  scale_size_manual(values = c(0.13, 0.29, 0.29)) +
  scale_colour_manual(values = c("#508CC8",
                                 "#FF4A4F",
                                 "#2D506E")) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    size = 0.4,
    color = "dark grey"
  ) +
  facet_grid(family ~ outcome, scales = "free", space = "free_y") +
  labs(x = "Beta, 95% CI (ml)", y = "") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    # remove y-axis info (as will be presented with brain volume)
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    strip.text.x =  element_text(size = 10.5,  face = "bold" , hjust = 0.1),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )

### Plot - white matter lesions ----
fig2a_wmhv <- fig2a_res %>%
  # filter for white matter lesions outcome
  filter(outcome == "C) White matter lesions") %>%
  ggplot(aes(y = pathogen)) +
  geom_pointrangeh(
    aes(
      x = beta,
      xmin = ci_lower,
      xmax = ci_upper,
      col = model,
      shape = threshold,
      size = threshold
    ),
    position = position_dodgev(height = 0.75)
  ) +
  scale_x_continuous(n.breaks = 6) +
  scale_y_discrete(limits = rev) +
  scale_shape_manual(values = c(4, 1, 19)) +
  scale_size_manual(values = c(0.13, 0.29, 0.29)) +
  scale_colour_manual(values = c("#508CC8",
                                 "#FF4A4F",
                                 "#2D506E")) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    size = 0.4,
    color = "dark grey"
  ) +
  facet_grid(family ~ outcome, scales = "free", space = "free_y") +
  labs(x = "% difference, 95% CI", y = "") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    # remove y-axis info (as will be presented with brain volume)
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    strip.text.x =  element_text(size = 10.5,  face = "bold" , hjust = 0.1),
    strip.placement = "outside",
    strip.background = element_blank(),
    # format facet text
    strip.text.y = element_text(
      face = "italic",
      size = 10,
      angle = 0,
      hjust = 0.1
    ),
    legend.position = "none"
  )

# Add lines for each pathogen family
fig2a_wmhv <- ggplotGrob(fig2a_wmhv)
lg <- linesGrob(
  x = unit(c(0, 0), "npc"),
  y = unit(c(0, 1), "npc"),
  gp = gpar(col = "black", lwd = 3)
)

for (k in grep("strip-r", fig2a_wmhv$layout$name)) {
  fig2a_wmhv$grobs[[k]]$grobs[[1]]$children[[1]] <- lg
}

grid.draw(fig2a_wmhv)

### Plot together & save as fig 2a ----
fig2a <- plot_grid(fig2a_bv,
                   fig2a_hv,
                   fig2a_wmhv,
                   nrow = 1,
                   rel_widths = c(1.4, 1, 1.7))

ggsave(
  plot = fig2a,
  filename = "figures/figure2a.png",
  height = 5.5,
  width = 9
)

## Figure 2b ----
### Format data ----
# Results data
fig2b_res <-
  # read in data
  read.csv("results/serostatus_main_MA.csv") %>%
  # filter for PBI results only
  filter(pathogen %in% c("PBI", "neuro_PBI")) %>%
  # tidy pbi name & convert to per 0.1 unit increase
  mutate(
    pathogen = case_when(as.character(pathogen) == "neuro_PBI" ~ "Neurotropic PBI \n(per 0.1 unit increase)",
                         as.character(pathogen) == "PBI" ~ "PBI \n(per 0.1 unit increase)"),
    # convert to per 0.1 unit increase
    beta = beta / 10,
    ci_lower = ci_lower / 10,
    ci_upper = ci_upper / 10
  ) %>%
  # convert to sympercents for white matter lesions
  mutate(
    beta = ifelse(outcome == "WMHV", beta * 100 , beta),
    ci_lower = ifelse(outcome == "WMHV", ci_lower * 100  , ci_lower),
    ci_upper = ifelse(outcome == "WMHV", ci_upper * 100  , ci_upper),
    # capitalise model for presentation purposes
    model = factor(str_to_sentence(model),
                   # set order for plotting
                   levels = c("Model 3", "Model 2", "Model 1")),
    # indicate significance threshold
    threshold = factor(
      case_when(
        padj < 0.05 ~ "padj<0.05",
        padj > 0.05 & p < 0.05 ~ "p<0.05",
        p > 0.05 ~ "p>0.05"
      ),
      levels = c("p>0.05",
                 "p<0.05",
                 "padj<0.05")
    ),
    # indicate family
    family = "Pathogen burden"
  ) %>%
  # select cols of interest
  select(pathogen,
         family,
         outcome,
         model,
         beta,
         p,
         ci_lower,
         ci_upper,
         threshold)

### Plot - brain volume ----
fig2b_bv <- fig2b_res %>%
  # filter for brain volume outcome
  filter(outcome == "brain volume") %>%
  # plot pathogen on y axis
  ggplot(aes(y = pathogen)) +
  # plot vals on x axis,
  geom_pointrangeh(
    aes(
      x = beta,
      xmin = ci_lower,
      xmax = ci_upper,
      # color by model
      col = model,
      # shape and size by threshold
      shape = threshold,
      size = threshold
    ),
    # plot models below each other
    position = position_dodgev(height = 0.75)
  ) +
  # specify number of breaks on x-axis
  scale_x_continuous(n.breaks = 6) +
  # reverse plotting of pathogens
  scale_y_discrete(limits = rev) +
  # set shape for plotting estimates (indicates threshold)
  scale_shape_manual(values = c(4, 1, 19)) +
  # set size for plotting estimate (indicates threshold)
  scale_size_manual(values = c(0.13, 0.29, 0.29)) +
  # set colors for plotting models
  scale_colour_manual(values = c("#508CC8",
                                 "#FF4A4F",
                                 "#2D506E")) +
  # include vertical dashed line at 0
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    size = 0.4,
    color = "dark grey"
  ) +
  # split into panels indicating pathogen family
  facet_grid(family ~ ., scales = "free", space = "free_y") +
  # amend axis labels
  labs(x = "Beta, 95% CI (ml)", y = "") +
  # set theme
  theme_classic() +
  theme(
    # set axis text size
    axis.text = element_text(size = 10),
    # amend space between facet panels
    panel.spacing = unit(0.8, "lines"),
    # remove y strip text & background
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    # remove legend
    legend.position = "none"
  )

### Plot - hippocampal volume ----
fig2b_hv <- fig2b_res %>%
  # filter for hippocampal volume outcome
  filter(outcome == "hippocampal volume") %>%
  ggplot(aes(y = pathogen)) +
  geom_pointrangeh(
    aes(
      x = beta,
      xmin = ci_lower,
      xmax = ci_upper,
      col = model,
      shape = threshold,
      size = threshold
    ),
    position = position_dodgev(height = 0.75)
  ) +
  scale_x_continuous(n.breaks = 6) +
  scale_y_discrete(limits = rev) +
  scale_shape_manual(values = c(4, 1, 19)) +
  scale_size_manual(values = c(0.13, 0.29, 0.29)) +
  scale_colour_manual(values = c("#508CC8",
                                 "#FF4A4F",
                                 "#2D506E")) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    size = 0.4,
    color = "dark grey"
  ) +
  facet_grid(family ~ ., scales = "free", space = "free_y") +
  labs(x = "Beta, 95% CI (ml)", y = "") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    # remove y-axis info (as will be presented with brain volume)
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0.81, "lines"),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )

### Plot - white matter lesions ----
fig2b_wmhv <- fig2b_res %>%
  # filter for white matter lesions outcome
  filter(outcome == "WMHV") %>%
  ggplot(aes(y = pathogen)) +
  geom_pointrangeh(
    aes(
      x = beta,
      xmin = ci_lower,
      xmax = ci_upper,
      col = model,
      shape = threshold,
      size = threshold
    ),
    position = position_dodgev(height = 0.75)
  ) +
  scale_x_continuous(n.breaks = 6) +
  scale_y_discrete(limits = rev) +
  scale_shape_manual(values = c(4, 1, 19)) +
  scale_size_manual(values = c(0.13, 0.29, 0.29)) +
  scale_colour_manual(values = c("#508CC8",
                                 "#FF4A4F",
                                 "#2D506E")) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    size = 0.4,
    color = "dark grey"
  ) +
  facet_grid(family ~ ., scales = "free", space = "free_y") +
  labs(x = "% difference, 95% CI", y = "") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    # remove y-axis info (as will be presented with brain volume)
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0.81, "lines"),
    strip.placement = "outside",
    strip.background = element_blank(),
    # format facet text
    strip.text.y = element_text(
      face = "italic",
      size = 10,
      angle = 0,
      hjust = 0.1
    ),
    legend.position = "none"
  )

# Add lines for each pathogen family
fig2b_wmhv <- ggplotGrob(fig2b_wmhv)
lg <- linesGrob(
  x = unit(c(0, 0), "npc"),
  y = unit(c(0, 1), "npc"),
  gp = gpar(col = "black", lwd = 3)
)

for (k in grep("strip-r", fig2b_wmhv$layout$name)) {
  fig2b_wmhv$grobs[[k]]$grobs[[1]]$children[[1]] <- lg
}

grid.draw(fig2b_wmhv)

### Plot together & save as fig 2b ----
fig2b <- plot_grid(
  fig2b_bv,
  fig2b_hv,
  fig2b_wmhv,
  nrow = 1,
  rel_widths = c(1.615, 1, 1.5)
)

## Figure 2 ----
ggsave(
  plot = plot_grid(
    fig2a,
    fig2b,
    nrow = 2,
    rel_heights = c(1, 0.205),
    align = "v",
    axis = "lr"
  ),
  filename = "figures/figure2.png",
  height = 8.25,
  width = 9.9
)

## End ----
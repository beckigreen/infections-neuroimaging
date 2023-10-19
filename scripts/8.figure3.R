# ============= #
#  8. Figure 3  #
# ============= #

## Set up ----
filepath <- ""
setwd(filepath)

## Load libraries ----
library(ggstance)
library(grid)
library(tidyverse)

## Format data ----
# List of antigens
antigens <- c(
  "HSV1 (1gG)",
  "HSV2 (2mgG)",
  "VZV (gEgI)",
  "EBV (VCA)*",
  "CMV (pp150)*",
  "HHV6A (1E1A)",
  "HHV6B (1E1B)",
  "HHV7 (U14)",
  "BK (VP1)",
  "JC (VP1)",
  "MCV (VP1)",
  "T.gondii (p22)*",
  "C.trachomatis (pGP3)",
  "H.pylori (CagA)"
)


# Results data
fig3_res <-
  # read in data
  read.csv("results/seroreac_main_MA.csv") %>%
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
    # format antigen names
    antigen = case_when(
      as.character(pathogen) == "HSV1_1gG_seroreac"  ~ "HSV1 (1gG)",
      as.character(pathogen) == "HSV2_2mgGunique_seroreac"  ~ "HSV2 (2mgG)",
      as.character(pathogen) == "VZV_gE_gI_seroreac"  ~ "VZV (gEgI)",
      as.character(pathogen) == "EBV_VCAp18_seroreac"  ~ "EBV (VCA)*",
      as.character(pathogen) == "CMV_pp150NTerm_seroreac"  ~ "CMV (pp150)*",
      as.character(pathogen) == "HHV6_1E1A_seroreac"  ~ "HHV6A (1E1A)",
      as.character(pathogen) == "HHV6_1E1B_seroreac"  ~ "HHV6B (1E1B)",
      as.character(pathogen) == "HHV7_U14_seroreac"  ~ "HHV7 (U14)",
      as.character(pathogen) == "BK_VP1_seroreac"  ~ "BK (VP1)",
      as.character(pathogen) == "JC_VP1_seroreac"  ~ "JC (VP1)",
      as.character(pathogen) == "MCV_VP1_seroreac"  ~ "MCV (VP1)",
      as.character(pathogen) == "Tg_p22_seroreac"  ~ "T.gondii (p22)*",
      as.character(pathogen) == "HP0547_1_seroreac"  ~ "H.pylori (CagA)",
      as.character(pathogen) == "Ct_pGP3_seroreac"  ~ "C.trachomatis (pGP3)"
    ),
    # format outcome name
    outcome = case_when(
      outcome == "brain volume" ~ "A) Brain volume",
      outcome == "hippocampal volume" ~ "B) Hippocampal volume",
      outcome == "WMHV" ~ "C) White matter lesions"
    )
  ) %>%
  mutate(
    # set order
    antigen = factor(antigen, levels = antigens),
    # indicate antigen family
    family = case_when(
      antigen %in% c(
        "HSV1 (1gG)",
        "HSV2 (2mgG)",
        "VZV (gEgI)",
        "EBV (VCA)*",
        "CMV (pp150)*",
        "HHV6A (1E1A)",
        "HHV6B (1E1B)",
        "HHV7 (U14)"
      ) ~ "Herpesviruses",
      antigen %in% c("BK (VP1)", "JC (VP1)", "MCV (VP1)") ~ "Human Polyomaviruses",
      antigen %in% c("T.gondii (p22)*", "C.trachomatis (pGP3)", "H.pylori (CagA)") ~ "Bacteria/Protozoa"
    )
  ) %>%
  # select cols of interest
  select(antigen,
         family,
         outcome,
         model,
         beta,
         p,
         ci_lower,
         ci_upper,
         threshold)

## Plot - brain volume ----
fig3_bv <- fig3_res %>%
  # filter for brain volume outcome
  filter(outcome == "A) Brain volume") %>%
  # plot antigen on y axis
  ggplot(aes(y = antigen)) +
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
  # reverse plotting of antigens
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
  # split into panels indicating antigen family
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

## Plot - hippocampal volume ----
fig3_hv <- fig3_res %>%
  # filter for hippocampal volume outcome
  filter(outcome == "B) Hippocampal volume") %>%
  ggplot(aes(y = antigen)) +
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

## Plot - white matter lesions ----
fig3_wmhv <- fig3_res %>%
  # filter for white matter lesions outcome
  filter(outcome == "C) White matter lesions") %>%
  ggplot(aes(y = antigen)) +
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

# Add lines for each antigen family
fig3_wmhv <- ggplotGrob(fig3_wmhv)
lg <- linesGrob(
  x = unit(c(0, 0), "npc"),
  y = unit(c(0, 1), "npc"),
  gp = gpar(col = "black", lwd = 3)
)

for (k in grep("strip-r", fig3_wmhv$layout$name)) {
  fig3_wmhv$grobs[[k]]$grobs[[1]]$children[[1]] <- lg
}

grid.draw(fig3_wmhv)

## Plot together & save as fig 3 ----
fig3 <- plot_grid(fig3_bv,
                  fig3_hv,
                  fig3_wmhv,
                  nrow = 1,
                  rel_widths = c(1.4, 1, 1.7))

ggsave(
  plot = fig3,
  filename = "figures/figure3.png",
  height = 5.3,
  width = 9
)

## End ----